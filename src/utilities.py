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
	"""
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
	
	for sample in projInfo.samples:
		binPath = binContigPath + '/' + sample
		if not os.path.exists(binPath):
			os.mdkir(binPath)
			
	contigIDs = {}
	for coreID in cores:
		for contigID in cores[coreID]:
			contigIDs[contigID] = coreID
			
	# create file handles
	ofhs = {}
	for sample in projInfo.samples:
		if sample not in ofhs:
			ofhs[sample] = {}
		for coreID in cores:
			binContigFile = binContigPath + '/' + sample + '/' + 
			if coreID not in ofhs[sample]:
				
	
	for sample in projInfo.samples:
			contigIDs[sample] = []
			assemblyFile = projInfo.getAssemblyFile(sample)
			afh = open(assemblyFile, 'r')
			for record in SeqIO.parse(afh, "fasta"):
				contigIDs[sample].append(record.id)
			contigIDs[sample] = set(contigIDs[sample])
			afh.close()
	"""
	sys.stdout.write('Code under construction\n')

def extractReadsForBins(projInfo, options):
	sys.stdout.write('Code under construction\n')