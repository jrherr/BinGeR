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
import cPickle
import pysam

def outputBins(projInfo, options):
	
	binContigPath = projInfo.out_dir + '/binContigs'
	
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
	for sample in projInfo.samples:
		if not options.quiet:
			sys.stdout.write('[%s]\r' % sample)
		ofhs = {}
		for coreID in cores:
			binContigFile = binContigPath + '/'+ coreID + '/' + sample + '.contigs.fa'
			ofh = open(binContigFile, 'w')
			ofhs[coreID] = ofh
			
		contigIDs[sample] = []
		assemblyFile = projInfo.getAssemblyFile(sample)
		afh = open(assemblyFile, 'r')
		for record in SeqIO.parse(afh, "fasta"):
			if record.id not in contigIDs:
				continue
			coreID = contigIDs[record.id]
			ofh = ofhs[coreID][sample]
			ofh.write('>%s\n%s\n'%(record.id, record.seq))
		afh.close()
	
		for coreID in ofhs:
			ofhs[coreID].close()
	
	sys.stdout.flush()
	if not options.quiet:
		sys.stdout.write('Done. Contigs stored at:\n %s\n' % binContigPath)
	
# End of outputBins

def extractReadsForBins(projInfo, options):
	binReadPath = projInfo.out_dir + '/binReads'
	
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
		sys.stdout.write('Now outputting contigs for each core...\n')
	
	# contigID lookup
	contigIDs = {}
	for coreID in cores:
		for contigID in cores[coreID]:
			contigIDs[contigID] = coreID			

	for sample in projInfo.samples:
		if not options.quiet:
			sys.stdout.write('[%s]\r' % sample)
		ofhs = {}
		for coreID in cores:
			binPEReadFile = binReadPath + '/'+ coreID + '/' + sample + '.PE.fa'
			binSEReadFile = binReadPath + '/'+ coreID + '/' + sample + '.SE.fa'
			ofh1 = open(binPEReadFile, 'w')
			ofh2 = open(binSEReadFile, 'w')
			ofhs[coreID] = (ofh1, ofh2)
		
		bamFile = projInfo.getBamFile(sample)
		samfh = pysam.Samfile(bamFile, 'rb')
		contigs = samfh.references()
		PEReadLookup = {}
		SEReadLookup = {}
		for contigID in contigs:
			if contigID not in contigIDs:
				continue
			coreID = contigIDs[contigID]
			readIDs = samfh.fetch(contigID)
			PEs, SEs = categorizeReads(readIDs)
			for x in PEs: PEReadLookup[x] = coreID
			for x in SEs: SEReadLookup[x] = coreID
		samfh.close()
		
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
				ofh = ofhs[coreID][sample][0]
				ofh.write('>%s\n%s\n' % (tag, seq))
			elif tag in SEReadLookup:
				coreID = SEReadLookup[tag]
				ofh = ofhs[coreID][sample][1]
				ofh.write('>%s\n%s\n' % (tag, seq))
		rfh.close()
	
		for sample in projInfo.samples:
			ofhs[coreID][0].close()
			ofhs[coreID][1].close()

	sys.stdout.flush()
	if not options.quiet:
		sys.stdout.write('Done. Contigs stored at:\n %s\n' % binContigPath)
	
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

# End of categorizeReads