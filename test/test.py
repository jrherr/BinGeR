#!/usr/bin/python
# -*- coding: utf-8 -*- 

__author__ = 'Chengwei Luo (luo.chengwei@gatech.edu)'
__version__ = '0.0.1'
__date__ = 'August 2013'

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

USAGE = \
"""Usage: %prog <required_parameters> [options]

BinGeR: in-situ Genome recovery tool for series metagenomes

BinGeR (Binner for Genome Recovery) is a platform for de novo genome recovery 
from series metagenomes. It integrates information from single copy gene content, 
contig coverage correlation, and contig oligo-nucleotide composition correlation, 
and contig alignments to bin genome fragments together. 
It is written in Python, therefore it should run on Windows, Mac OS, and Unix/Linux. 

Add --help to see a full list of required and optional
arguments to run BinGeR.

Additional information can also be found at:
https://github.com/luo-chengwei/BinGeR/wiki

If you use BinGeR in your work, please cite it as:
<BinGeR citation here>

Copyright: Chengwei Luo, Konstantinidis Lab, Georgia Institute of Technology, 2013
"""



import sys
import re
import glob
import os
from optparse import OptionParser, OptionGroup
import networkx as nx
import cPickle
import tarfile
import shutil
from subprocess import PIPE, Popen

import taxonomy
import phylo

class ProjectInfo:
	def __init__(self):
		self.samples = []
		self.reads_dir = None
		self.bams_dir = None
		self.coverage_dir = None
		self.assemblies_dir = None
		self.zscore_dir = None
		self.hmmscan_dir = None
		self.out_dir = None
		self.HMM = {}
		self.DBs = {}
		
	def initPaths(self, options):
		if os.path.exists(options.assemblies_dir):
			self.assemblies_dir = os.path.abspath(options.assemblies_dir)
		else:
			sys.stderr.write("FATAL: cannot find assemblies directory:%s\n"%options.reads_dir)
			exit(0)
		
		if os.path.exists(options.out_dir):
			self.out_dir = os.path.abspath(options.out_dir)
		else:
			try:
				os.mkdir(options.out_dir)
			except OSError:
				sys.stderr.write("FATAL: cannot create output directory %s, please make sure you have proper privilege!\n"%options.out_dir)
				exit(0)
			
	def getAssemblyFile(self, sample):
		assemblyFiles = glob.glob(self.assemblies_dir+'/'+sample+'[.]*[.]fa')
		if len(assemblyFiles) == 0:
			assemblyFiles = glob.glob(self.assemblies_dir+'/'+sample+'[.]fa')
		if len(assemblyFiles) == 0:
			sys.stderr.write("FATAL: fail to locate assembly for sample:%s\n"%sample)
			exit(0)
		elif len(assemblyFiles) > 1:
			sys.stderr.write("FATAL: ambiguity in locating assembly for sample:%s\n"%sample)
			sys.stderr.write("\tBinGeR found:\n")
			for readFile in assemblyFiles:
				sys.stderr.write("\t%s\n"%readFile)
			exit(0)
		else:
			return assemblyFiles[0]
			
	def initProject(self, options):
		sys.stdout.write("Fetching information for this BinGeR project...\n")

		# Reads the sample list and check if every required file is in place
		try:
			sfh = open(options.sample_list, 'r')
		except:
			sys.stderr.write("FATAL: Cannot open supplied sample list: %s\n"%options.sample_list)
			exit(0)
			
		while 1:
			sample = sfh.readline().rstrip('\n')
			if not sample:
				break
			sample = sample.replace(' ', '', sample.count(' '))
			if sample != '':
				self.samples.append(sample)
		sfh.close()
		
		self.initPaths(options)
		
		sys.stdout.write("Checking if every file in place...\n")
		for sample in self.samples:
			assemblyFile = self.getAssemblyFile(sample)
			
		sys.stdout.write("Done fetching information!\n")
		
	
		# check if db/ has every file needed.
		sys.stdout.write("Checking db files...\n")
		db_dir = os.path.abspath(sys.argv[0]).replace('test.py', 'db/')
		dbs = ['HMM.txt', 'ncbiNodes.lib','ncbiSciNames.lib', 
				'singleCopy.prot.tar.gz', 'singleCopy.nuc.tar.gz']
		dbFiles = [db_dir + filename for filename in dbs]
		for dbFile in dbFiles:
			if not os.path.exists(dbFile):
				sys.stderr.write('FATAL: Cannot locate db file: %s\n' % dbFile)
				exit(0)
		
		self.DBs['HMM'] = dbFiles[0]
		self.DBs['ncbiNodes'] = dbFiles[1]
		self.DBs['ncbiSciNames'] = dbFiles[2]
		self.DBs['nuc'] = dbFiles[3]
		self.DBs['prot'] = dbFiles[4]
		sys.stdout.write("Done. Everything looks fine.\n")
		
		
		# check db/HMM.txt and load it to self.HMM
		sys.stdout.write("Loading HMMFams...\n")
		hmmFile = self.DBs['HMM']
		if not os.path.exists(hmmFile):
			sys.stderr.write("FATAL: cannot locate db/HMM.txt, abort BinGeR.\n")
			exit(0)
		
		hmmfh = open(hmmFile, 'rb')
		self.HMM = cPickle.load(hmmfh)
		hmmfh.close()
		sys.stdout.write('HMMFams loaded.\n')
		
# end of class ProjectInfo


def main(argv = sys.argv[1:]):

	parser = OptionParser(usage = USAGE, version="Version: " + __version__)
	
	# Required arguments
	requiredOptions = OptionGroup(parser, "Required options",
								"These options are required to run BinGeR, and may be supplied in any order.")
	
	requiredOptions.add_option("-l", "--sample_list", type = "string", metavar = "FILE",
							help = "Text file containing all sample names, one per line")

	requiredOptions.add_option("-o", "--out_dir", type = "string", metavar = "OUTDIR",
							help = "Working directory where the results and intermediate files will be stored at")
	
	requiredOptions.add_option("-i", "--core", type = "string", metavar = "FILE",
							help = "Input init core")

	parser.add_option_group(requiredOptions)

	# Optional arguments that need to be supplied if not the same as default
	optOptions = OptionGroup(parser, "Optional parameters",
						"There options are optional, and may be supplied in any order.")

	optOptions.add_option("-a", "--assemblies_dir", type = "string", default = "Assemblies", metavar = "DIR",
							help = "Directory where assemblies in fasta format are, naming follows \"sample.*.fa\" convention. [Default: ./Assemblies]")

	optOptions.add_option("--blat", type = "string", default = "blat",
							help = "Path to blat, specify if not in env.")

	parser.add_option_group(optOptions)
	
	# runtime settings that could affect the file saving and message printing
	runtimeSettings = OptionGroup(parser, "Runtime settings",
						"There options are optional, and may be supplied in any order.")
						
	runtimeSettings.add_option("-q", "--quiet", default = False, action = "store_true",
								help = "Suppress printing detailed runtime information, only important messages will show [default: False].")

	runtimeSettings.add_option("--no_intermediates", action="store_false", 
								dest = "save_intermediates", default = True,
								help = "Do no save intermediate files during runtime [default: True (save intermediates)].")

	parser.add_option_group(runtimeSettings)
	
	(options, args) = parser.parse_args(argv)
	
	# test if blat exists
	blatTest = Popen(options.blat, shell=True, stdout=PIPE).stdout.read()
	if blatTest == None or len(blatTest) == 0:
		sys.stderr.write("FATAL: blat not found in path!")
		exit(0)
	
	# check sanity of the files in required directories
	projInfo = ProjectInfo()
	projInfo.initProject(options)

	sys.stdout.write('Unpickling core...\n')
	ifh = open(options.core, 'rb')
	G = cPickle.load(ifh)
	ifh.close()
	sys.stdout.write('Done.\n')
	
	sys.stdout.write('Loading tTree...\n')
	tTree = taxonomy.TaxonTree()
	tTree.loadTreeFromNodeLib(projInfo.DBs['ncbiNodes'], projInfo.DBs['ncbiSciNames'])
	sys.stdout.write('Done.\n')
	
	# extract db files to out_dir
	try:
		nucTar = tarfile.open(projInfo.DBs['nuc'], 'r')
		protTar = tarfile.open(projInfo.DBs['prot'], 'r')
	except:
		sys.stderr.write('FATAL: failure in opening tarfile.\n')
		exit(0)
	
	nucTar.extractall(path = projInfo.out_dir)
	protTar.extractall(path = projInfo.out_dir)
	
	"""
	try:
		nucTar.extractall(path = projInfo.out_dir)
	except:		
		sys.stderr.write('FATAL: failure in extracting single copy genes nucleotide sequences.\n')
		exit(0)
			
	try:
		protTar.extractall(path = projInfo.out_dir)
	except:		
		sys.stderr.write('FATAL: failure in extracting single copy genes protein sequences.\n')
		exit(0)
	"""
	
	# create merged files
	catNuc = open(projInfo.out_dir + '/genes.nuc.clustered/allGenes.fasta', 'w')
	catProt = open(projInfo.out_dir + '/genes.prot.clustered/allGenes.fasta', 'w')
	
	for file in glob.glob(projInfo.out_dir + '/genes.nuc.clustered/*.fa'):
		shutil.copyfileobj(open(file, 'rb'), catNuc)
	catNuc.close()
		
	for file in glob.glob(projInfo.out_dir + '/genes.prot.clustered/*.fa'):
		shutil.copyfileobj(open(file, 'rb'), catProt)
	catProt.close()

	P = phylo.nodePhylo(G, tTree, projInfo, options)
	
	# clean up the directory
#	shutil.rmtree(projInfo.out_dir + '/genes.nuc.clustered/')
#	shutil.rmtree(projInfo.out_dir + '/genes.prot.clustered/')
	"""
	try:
		shutil.rmtree(projInfo.out_dir + '/genes.nuc.clustered/')
	except:
		sys.stderr.write('WARNING: Failure in removing intermediate nucleotide sequences directory.\n')
			
	try:
		shutil.rmtree(projInfo.out_dir + '/genes.prot.clustered/')
	except:
		sys.stderr.write('WARNING: Failure in removing intermediate protein sequences directory.\n')	
	"""

if __name__ == '__main__':
	main()