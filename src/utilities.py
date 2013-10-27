#!/usr/bin/python
# -*- coding: utf-8 -*- 


"""
BinGeR (Binner for Genome Recovery): 
	in-situ Genome recovery tool for series metagenomes

Copyright(c) 2013 Chengwei Luo (luo.chengwei@gatech.edu)

	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

https://github.com/luo-chengwei/BinGeR

for help, type:
python BinGeR.py --help
"""

import sys
import os
import glob
import cPickle
import shutil
import resource
from os.path import realpath
from operator import itemgetter

import pysam
from Bio import SeqIO

NOFILE_LIMIT, NOVFILE_LIMIT = resource.getrlimit(resource.RLIMIT_NOFILE)

def outputBins(projInfo, options):
	binContigPath = projInfo.out_dir + '/binContigs'
	if os.path.exists(binContigPath):
		if len(glob.glob(binContigPath + '/*')) > 0:
			return
		else:
			pass		
	else:
		os.mkdir(binContigPath)

	finalCoresPickle = projInfo.out_dir + '/finalCores.cpickle'
	if not os.path.exists(finalCoresPickle):
		sys.stderr.write('FATAL: failure in locating the final cores serialized results.\n')
		exit(0)
	
	try:
		pfh = open(finalCoresPickle, 'rb')
		cores = cPickle.load(pfh)
		pfh.close()
	except:
		sys.stderr.write('FATAL: failure in unpickling the final cores.\n')
		exit(0)
	
	for coreID in cores:
		binPath = binContigPath + '/' + coreID
		if not os.path.exists(binPath):
			os.mkdir(binPath)
			
	contigIDs = {}
	for coreID in cores:
		for contigID in cores[coreID]:
			contigIDs[contigID] = coreID
	
	if not options.quiet:
		sys.stdout.write('Now outputting contigs for each core...\n')
	
	# dynamic file handle dict
	numFileHandles = len(projInfo.samples) * len(cores)
	ofhs = []
	for i in range(numFileHandles):
		ofhs.append(None)
	
	coreIDs = cores.keys()
	activeFileHandles = 0
	for i, sample in enumerate(projInfo.samples):
		for j, coreID in enumerate(coreIDs):
			ofhIndex = j + (i * len(projInfo.samples))
			binContigFile = binContigPath + '/'+ coreID + '/' + sample + '.contigs.fa'
			if activeFileHandles < NOFILE_LIMIT - len(projInfo.samples) - 4:
				ofhs[ofhIndex] = open(binContigFile, 'a')
				activeFileHandles += 1
			else:
				ofhs[ofhIndex] = None
	
	for i, sample in enumerate(projInfo.samples):
		if not options.quiet:
			sys.stdout.write('[%s]\n' % sample)
		
		contigIDs[sample] = []
		assemblyFile = projInfo.getAssemblyFile(sample)
		afh = open(assemblyFile, 'r')
		unclassifiedContigs = binContigPath + '/' + sample + '.unclassified.contigs.fa'
		ufh = open(unclassifiedContigs, 'w')
		
		for record in SeqIO.parse(afh, "fasta"):
			if record.id not in contigIDs:
				ufh.write('>%s\n%s\n'%(record.id, record.seq))
				continue
			coreID = contigIDs[record.id]
			j = coreIDs.index(coreID)
			
			# dynamically control filehandles, so the total # of active 
			# handles is smaller than the resource limit
			
			ofhIndex = j + (i * len(projInfo.samples))
			if ofhs[ofhIndex] == None:
				binContigFile = binContigPath + '/'+ coreID + '/' + sample + '.contigs.fa'
				# choose an active one to close
				fhSuc = 0
				for index, ofh in enumerate(ofhs[ofhIndex+1:]):
					if ofh != None:
						ofhs[index + ofhIndex + 1].close()
						ofhs[index + ofhIndex + 1] = None
						fhSuc = 1
						break
				if fhSuc == 0:
					for index, ofh in enumerate(ofhs[:ofhIndex]):
						if ofh != None:
							ofhs[index].close()
							ofhs[index] = None
							fhSuc = 1
							break
				# and then open the one we need
				ofhs[ofhIndex] = open(binContigFile, 'a')
			else:
				pass
			
			ofhs[ofhIndex].write('>%s\n%s\n'%(record.id, record.seq))
		
		ufh.close()
		afh.close()
	
	# close up all filehandles
	for ofh in ofhs:
		if ofh == None:
			continue	
		ofh.close()
	
	sys.stdout.flush()
	if not options.quiet:
		sys.stdout.write('Done. Contigs stored at:\n %s\n' % binContigPath)
	
# End of outputBins

def extractReadsForBins(projInfo, options):
	binReadPath = projInfo.out_dir + '/binReads'
	if os.path.exists(binReadPath):
		if len(glob.glob(binReadPath + '/*')) > 0:		
			return
		else:
			pass
	else:
		os.mkdir(binReadPath)
		
	finalCoresPickle = projInfo.out_dir + '/finalCores.cpickle'
	if not os.path.exists(finalCoresPickle):
		sys.stderr.write('FATAL: failure in locating the final cores serialized results.\n')
		exit(0)
	
	try:
		pfh = open(finalCoresPickle, 'rb')
		cores = cPickle.load(pfh)
		pfh.close()
	except:
		sys.stderr.write('FATAL: failure in unpickling the final cores.\n')
		exit(0)
	
	for coreID in cores:
		binPath = binReadPath + '/' + coreID
		if not os.path.exists(binPath):
			os.mkdir(binPath)
	
	if not options.quiet:
		sys.stdout.write('Now outputting reads for each core...\n')
	
	# dynamic file handle dict
	numFileHandles = 2 * len(projInfo.samples) * len(cores)
	ofhs = []
	for i in range(numFileHandles):
		ofhs.append(None)
	
	coreIDs = cores.keys()
	activeFileHandles = 0
	for i, sample in enumerate(projInfo.samples):
		for j, coreID in enumerate(coreIDs):
			
			ofhIndex1 = 2 * (j + (i * len(projInfo.samples)))
			ofhIndex2 = ofhIndex1 + 1
			binPEReadFile = binReadPath + '/'+ coreID + '/' + sample + '.PE.fa'
			binSEReadFile = binReadPath + '/'+ coreID + '/' + sample + '.SE.fa'
			
			if activeFileHandles < NOFILE_LIMIT -6:
				ofhs[ofhIndex1] = open(binPEReadFile, 'a')
				ofhs[ofhIndex2] = open(binSEReadFile, 'a')
				activeFileHandles += 2
			else:
				ofhs[ofhIndex1] = None
				ofhs[ofhIndex2] = None
	
	# contigID lookup
	contigIDs = {}
	for coreID in cores:
		for contigID in cores[coreID]:
			contigIDs[contigID] = coreID			

	for i, sample in enumerate(projInfo.samples):
		if not options.quiet:
			sys.stdout.write('Working on sample: %s\r' % sample)
			
		bamFile = projInfo.getBamFile(sample)
		
		try:
			samfh = pysam.Samfile(bamFile, 'rb')
		except IOError:
			sys.stderr.write('Failure in opening:\n\t%s\n' % bamFile)
			continue
		
		contigs = samfh.references
		PEReadLookup = {}
		SEReadLookup = {}
		for contigID in contigs:
			if contigID not in contigIDs:
				continue
			coreID = contigIDs[contigID]
			readIDs = []
			for read in samfh.fetch(contigID):
				if read.is_secondary:
					continue
				readIDs.append(read.qname)
			
			PEs, SEs = categorizeReads(readIDs)
			for x in PEs: PEReadLookup[x] = coreID
			for x in SEs: SEReadLookup[x] = coreID
			print PEs[:4]
		
		samfh.close()
		print 'Finsihed reading the samfile'
		
		readFile = projInfo.getReadFile(sample)
		rfh = open(readFile, 'r')
		while 1:
			tag = rfh.readline().rstrip('\n')
			if not tag:
				break
			tag = tag.replace('>', '')
			seq = rfh.readline().rstrip('\n')
			
			if tag not in PEReadLookup and tag not in SEReadLookup:
				continue
			elif tag in PEReadLookup:
				coreID = PEReadLookup[tag]
				j = coreIDs.index(coreID)
				ofhIndex = 2 * (j + (i * len(projInfo.samples)))
			elif tag in SEReadLookup:
				coreID = SEReadLookup[tag]
				j = coreIDs.index(coreID)
				ofhIndex = 1+ 2 * (j + (i * len(projInfo.samples)))
			
			if ofhs[ofhIndex] == None:
				if ofhIndex % 2 == 0:
					readFile = binReadPath + '/'+ coreID + '/' + sample + '.PE.fa'
				else:
					readFile = binReadPath + '/'+ coreID + '/' + sample + '.SE.fa'
				
				# choose an active one to close
				fhSuc = 0
				for index in range(ofhIndex+2, len(ofhs), 2):
					if ofhs[index] != None:
						ofhs[index].close()
						ofhs[index] = None
						fhSuc = 1
						break
				if fhSuc == 0:
					if ofhIndex % 2 == 0:
						startIndex = 0
					else:
						startIndex = 1
					for index in range(startIndex, ofhIndex, 2):
						if ofhs[index] != None:
							ofhs[index].close()
							ofhs[index] = None
							fhSuc = 1
							break
				# and then open the one we need
				ofhs[ofhIndex] = open(readFile, 'a')
				
			else:
				ofhs[ofhIndex].write('>%s\n%s\n' % (tag, seq))
		
		rfh.close()
		
		if not options.quiet:
			sys.stdout.flush()
			
	# close up all filehandles
	for ofh in ofhs:
		if ofh == None:
			continue	
		ofh.close()
		
	sys.stdout.flush()
	
	if not options.quiet:
		sys.stdout.write('Done. Reads stored at:\n %s\n' % binReadPath)
	
# End of extractReadsForBins

def categorizeReads(readIDs):
	PEs = []
	SEs = []
	occurrences = {}
	for readID in readIDs:
		if readID[:-1] not in occurrences:
			occurrences[readID[:-1]] = []
		occurrences[readID[:-1]].append(readID)
	
	for readID in occurrences:
		if len(occurrences[readID]) == 1:
			SEs += occurrences[readID]
		else:
			PEs += occurrences[readID]
	return PEs, SEs
# End of categorizeReads

def cleanup(projInfo):
	for pickle in glob.glob(projInfo.out_dir + '/*.cpickle'):
		os.remove(pickle)
	shutil.rmtree(projInfo.out_dir + '/blat')
	shutil.rmtree(projInfo.out_dir + '/initCores')
	shutil.rmtree(projInfo.out_dir + '/subgraphs')

# End fo cleanup