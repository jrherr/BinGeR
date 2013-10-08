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
import os
import re
import glob
from optparse import OptionParser, OptionGroup
from subprocess import PIPE, Popen
import cPickle

from time import ctime, time
from datetime import timedelta

import contigSpace as cSpace
import utilities
import phylo

################################## CLASSES ################################
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
		if os.path.exists(options.reads_dir):
			self.reads_dir = os.path.abspath(options.reads_dir)
		else:
			sys.stderr.write("FATAL: cannot find reads directory:%s\n"%options.reads_dir)
			exit(0)
		
		if os.path.exists(options.bams_dir):
			self.bams_dir = os.path.abspath(options.bams_dir)
		else:
			sys.stderr.write("FATAL: cannot find bam files directory:%s\n"%options.bams_dir)
			exit(0)
		
		if os.path.exists(options.coverage_dir):
			self.coverage_dir = os.path.abspath(options.coverage_dir)
		else:
			sys.stderr.write("FATAL: cannot find coverage files directory:%s\n"%options.coverage_dir)
			exit(0)
		
		if os.path.exists(options.assemblies_dir):
			self.assemblies_dir = os.path.abspath(options.assemblies_dir)
		else:
			sys.stderr.write("FATAL: cannot find assemblies directory:%s\n"%options.reads_dir)
			exit(0)
		
		if os.path.exists(options.hmmscan_dir):
			self.hmmscan_dir = os.path.abspath(options.hmmscan_dir)
		else:
			sys.stderr.write("FATAL: cannot find hmmscan files directory:%s\n"%options.hmmscan_dir)
			exit(0)
		
		if os.path.exists(options.zscore_dir):
			self.zscore_dir = os.path.abspath(options.zscore_dir)
		else:
			sys.stderr.write("FATAL: cannot find zscore files directory:%s\n"%options.zscore_dir)
			exit(0)
		
		if os.path.exists(options.out_dir):
			self.out_dir = os.path.abspath(options.out_dir)
		else:
			try:
				os.mkdir(options.out_dir)
			except OSError:
				sys.stderr.write("FATAL: cannot create output directory %s, please make sure you have proper privilege!\n"%options.out_dir)
				exit(0)
		
	def printSamples(self):
		sys.stdout.write("Samples used in this BinGeR run:\n")
		for sample in self.samples:
			sys.stdout.write("\t%s\n"%sample)
	
	def getBamFile(self, sample):
		bamFiles = glob.glob(self.bams_dir+'/'+sample+'[.]*[.]bam')
		if len(bamFiles) == 0:
			sys.stderr.write("FATAL: fail to locate bam file for sample:%s"%sample)
			exit(0)
		elif len(bamFiles) > 1:
			sys.stderr.write("FATAL: ambiguity in locating bam file for sample:%s"%sample)
			sys.stderr.write("\tBinGeR found:\n")
			for bamFile in bamFiles:
				sys.stderr.write("\t%s\n"%bamFile)
			exit(0)
		else:
			return bamFiles[0]
			
	def getReadFile(self, sample):
		readFiles = glob.glob(self.reads_dir+'/'+sample+'[.]*[.]fa')
		if len(readFiles) == 0:
			readFiles = glob.glob(self.reads_dir+'/'+sample+'[.]fa')
		if len(readFiles) == 0:
			sys.stderr.write("FATAL: fail to locate read file for sample:%s\n"%sample)
			exit(0)
		elif len(readFiles) > 1:
			sys.stderr.write("FATAL: ambiguity in locating reads for sample:%s\n"%sample)
			sys.stderr.write("\tBinGeR found:\n")
			for readFile in readFiles:
				sys.stderr.write("\t%s\n"%readFile)
			exit(0)
		else:
			return readFiles[0]
	
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
			
	def getCoverageFile(self, sampleA, sampleB):
		coverageFiles = glob.glob(self.coverage_dir+'/'+sampleA+'[.]vs[.]'+sampleB+'[.]*[.]coverage')
		if len(coverageFiles) == 0:
			coverageFiles = glob.glob(self.coverage_dir+'/'+sampleA+'[.]vs[.]'+sampleB+'[.]coverage')
		if len(coverageFiles) == 0:
			sys.stderr.write("FATAL: fail to locate coverage file for %s vs %s\n"%sampleA, sampleB)
			exit(0)
		elif len(coverageFiles) > 1:
			sys.stderr.write("FATAL: ambiguity in locating coverage file for %s vs %s\n"%sampleA, sampleB)
			sys.stderr.write("\tBinGeR found:\n")
			for coverageFile in coverageFiles:
				sys.stderr.write("\t%s\n"%coverageFile)
			exit(0)
		else:
			return coverageFiles[0]
			
	def getZScoreFile(self, sample):
		ZScoreFiles = glob.glob(self.zscore_dir+'/'+sample+'[.]*[.]ZScore')
		if len(ZScoreFiles) == 0:
			ZScoreFiles = glob.glob(self.zscore_dir+'/'+sample+'[.]ZScore')
		if len(ZScoreFiles) == 0:
			sys.stderr.write("FATAL: fail to locate z-score file for sample:%s\n"%sample)
			exit(0)
		elif len(ZScoreFiles) > 1:
			sys.stderr.write("FATAL: ambiguity in locating ZScore file for sample:%s\n"%sample)
			sys.stderr.write("\tBinGeR found:\n")
			for ZScoreFile in ZScoreFiles:
				sys.stderr.write("\t%s\n"%ZScoreFile)
			exit(0)
		else:
			return ZScoreFiles[0]
			
	def getHMMScanFile(self, sample):
		HMMScanFiles = glob.glob(self.hmmscan_dir+'/'+sample+'[.]*[.]hmmscan')
		if len(HMMScanFiles) == 0:
			HMMScanFiles = glob.glob(self.hmmscan_dir+'/'+sample+'[.]hmmscan')
		if len(HMMScanFiles) == 0:
			sys.stderr.write("FATAL: fail to locate hmmscan file for sample:%s\n"%sample)
			exit(0)
		elif len(HMMScanFiles) > 1:
			sys.stderr.write("FATAL: ambiguity in locating hmmscan file for sample:%s\n"%sample)
			sys.stderr.write("\tBinGeR found:\n")
			for HMMScanFile in HMMScanFiles:
				sys.stderr.write("\t%s\n"%HMMScanFile)
			exit(0)
		else:
			return HMMScanFiles[0]
			
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
			readFile = self.getReadFile(sample)
			bamFile = self.getBamFile(sample)
			assemblyFile = self.getAssemblyFile(sample)
			zscoreFile = self.getZScoreFile(sample)
			hmmscanFile = self.getHMMScanFile(sample)
			
		for sampleA in self.samples:
			for sampleB in self.samples:
				coverageFile = self.getCoverageFile(sampleA, sampleB)
				
		sys.stdout.write("Done fetching information!\n")
		
	
		# check if db/ has every file needed.
		sys.stdout.write("Checking db files...\n")
		db_dir = os.path.abspath(sys.argv[0]).replace('bin/BinGeR.py', 'db/')
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

################################### MAIN ######################################
def main(argv = sys.argv[1:]):

	parser = OptionParser(usage = USAGE, version="Version: " + __version__)
	
	# Required arguments
	requiredOptions = OptionGroup(parser, "Required options",
								"These options are required to run BinGeR, and may be supplied in any order.")
	
	requiredOptions.add_option("-l", "--sample_list", type = "string", metavar = "FILE",
							help = "Text file containing all sample names, one per line")

	requiredOptions.add_option("-o", "--out_dir", type = "string", metavar = "OUTDIR",
							help = "Working directory where the results and intermediate files will be stored at")

	parser.add_option_group(requiredOptions)

	# Optional arguments that need to be supplied if not the same as default
	optOptions = OptionGroup(parser, "Optional parameters",
						"There options are optional, and may be supplied in any order.")

	optOptions.add_option("-r", "--reads_dir", type = "string", default = "Reads/", metavar = "DIR",
							help = "Directory where coupled reads (interleaved) in fasta format are, the naming should follow \"sample.*.fa\" convention. [Default: ./Reads]")

	optOptions.add_option("-b", "--bams_dir", type = "string", default = "Bams", metavar = "DIR",
							help = "Directory where sorted bam files (reads versus assembly, same sample) are, the naming should follow \"sample.*.bam\" convention. [Default: ./Bams]")

	optOptions.add_option("-c", "--coverage_dir", type = "string", default = "Coverage", metavar = "DIR",
							help = "Directory where coverage files are, naming follows \"sampleA.vs.sampleB.*.coverage\" convention. [Default: ./Coverage]")

	optOptions.add_option("-a", "--assemblies_dir", type = "string", default = "Assemblies", metavar = "DIR",
							help = "Directory where assemblies in fasta format are, naming follows \"sample.*.fa\" convention. [Default: ./Assemblies]")

	optOptions.add_option("-z", "--zscore_dir", type = "string", default = "ZScores", metavar = "DIR",
							help = "Directory where oligo-nt z-score files are, naming follows \"sample.*.ZScore\" convention. [Default: ./ZScore]")

	optOptions.add_option("-s", "--hmmscan_dir", type = "string", default = "HMMScan", metavar = "DIR",
							help = "Directory where hmmscan files are, naming follows \"sample.*.hmmscan\" convention. [Default: ./HMMScan]")

	optOptions.add_option("-t", "--num_proc", type = "int", default = 1, metavar = 'INT',
							help = "Number of processor for BinGeR to use [default: 1].")
						
	optOptions.add_option("--blat", type = "string", default = "blat",
							help = "Path to blat, specify if not in env.")

	parser.add_option_group(optOptions)
	
	
	# Binning parameters that could fine tune the process
	clusteringOptions = OptionGroup(parser, "Binning parameters",
						"There options are optional, and may be supplied in any order.")
	
	clusteringOptions.add_option("-m", "--min_core", type = "int", default = 1e5, metavar = 'INT',
							help = "Minimum size to consider as bin core [default: 1e5].")
	
	clusteringOptions.add_option("-u", "--cov_clustering_min_length", dest = "minCovLength",
							type = "int", default = 1500, metavar = 'INT',
							help = "Minimum contig length to be considered in coverage clustering [default: 1500].")
	
	clusteringOptions.add_option("--min_cov_corrcoef", dest = "minCovCorrceof",
							type = "float", default = 0.95, metavar = 'FLOAT',
							help = "Minimum correlation coefficient cutoff for form a link between contigs using coverage profiles [default: 0.95].")
	
	clusteringOptions.add_option("--min_zscore_corrcoef", dest = "minZScoreCorrceof",
							type = "float", default = 0.9, metavar = 'FLOAT',
							help = "Minimum correlation coefficient cutoff for form a link between contigs using tri-/tetra-nt frequency Z-Score [default: 0.90].")
						
	clusteringOptions.add_option("-x", "--zscore_clustering_min_length", dest = "minZLength",
							type = "int", default = 3000, metavar = 'INT',
							help = "Minimum contig length to be considered in Z-score clustering [default: 2000].")
	
	
	clusteringOptions.add_option("--cpr_alpha", type = "float", default = 0.99, metavar = 'FLOAT',
							help = "The dampening factor, alpha, in community personalized PageRank [default: 0.99].")
	
	clusteringOptions.add_option("--cpr_tol", type = "float", default = 0.01, metavar = 'FLOAT',
							help = "The toll in community personalized PageRank [default: 0.01].")
							
	
	parser.add_option_group(clusteringOptions)

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
	
	if options.sample_list is None:
		parser.error("A list of samples and a working directory are required!")
		exit(0)
		
	if options.out_dir is None:
		parser.error("An output directory is required to supply!")
		exit(0)
		
	if options.num_proc < 1:
		parser.error("Number of processors must be integer >= 1, you supplied %i" % options.num_proc)
		exit(0)
		
	if options.min_core < 1e4 or options.min_core > 1e6:
		parser.error("Size of minimum bin size must be in range [1e4, 1e6]bp, you supplied %i" % options.min_core)
		exit(0)
	
	
	if options.cpr_alpha <= 0.5 or options.cpr_alpha >= 1:
		parser.error("Community PageRank Alpha must be in range (0.5, 1), you supplied %.3f" % options.cpr_alpha)
		exit(0)
	
	if options.cpr_tol <= 0 or options.cpr_tol >= 0.2:
		parser.error("Community PageRank Alpha must be in range (0, 0.2), you supplied %.3f" % options.cpr_tol)
		exit(0)

	total_start_time = time()
	sys.stdout.write("BinGeR started at %s\n"%(ctime()))
	sys.stdout.flush()

	# test if blat exists
	blatTest = Popen(options.blat, shell=True, stdout=PIPE).stdout.read()
	if blatTest == None or len(blatTest) == 0:
		sys.stderr.write("FATAL: blat not found in path!")
		exit(0)
	
	# check sanity of the files in required directories
	projInfo = ProjectInfo()
	projInfo.initProject(options)
	if not options.quiet:
		projInfo.printSamples()
	
	# build networkx graph for Project
	sys.stdout.write('\nInitializing contig space...\n')
	G = cSpace.ContigSpace(projInfo.samples)
	
	initCoresPath = projInfo.out_dir + '/initCores'
	refinedCoresPath = projInfo.out_dir + '/refinedCores'
	finalCoresPath = projInfo.out_dir + '/finalCores'
	if os.path.exists(finalCoresPath):
		if len(glob.glob(finalCoresPath + '/*.cpickle')) > 0:
			pass
	elif os.path.exists(refinedCoresPath):
		if len(glob.glob(refinedCoresPath + '/*.cpickle')) > 0:
			G.recruitContigs(projInfo, options)
	elif os.path.exists(initCoresPath):
		if len(glob.glob(initCoresPath + '/*.cpickle')) > 0:
			G.refineCores(projInfo, options)
			G.recruitContigs(projInfo, options)
	else:
		G.initSubgraphs(projInfo, options)
		G.forgeCores(projInfo, options)
		G.refineCores(projInfo, options)
		G.recruitContigs(projInfo, options)
	
	# output bins and the evaluation, extract reads of bins for downstream analysis.
	binContigPath = projInfo.out_dir + '/binContigs'
	if os.path.exists(binContigPath):
		if len(glob.glob(binContigPath + '/*.contigs.fa')) > 0:
			pass
		else:
			sys.stdout.write('Now outputting bins and statistics...\n')
			utilities.outputBins(projInfo, options)
			sys.stdout.write('Done.\n')
	else:
		sys.stdout.write('Now outputting bins and statistics...\n')
		utilities.outputBins(projInfo, options)
		sys.stdout.write('Done.\n')
	
	# get all the reads for the bins
	binReadPath = projInfo.out_dir + '/binReads'
	if os.path.exists(binReadPath):
		if len(glob.glob(binReadPath + '/*.reads.fa')) > 0:
			pass
		else:
			sys.stdout.write('Now extract reads for bins...\n')
			utilities.extractReadsForBins(G, projInfo, options)
			sys.stdout.write('Done.\n')
	else:
		sys.stdout.write('Now extract reads for bins...\n')
		utilities.extractReadsForBins(projInfo, options)
		sys.stdout.write('Done.\n')
	total_finish_time = time()
	sys.stdout.write("BinGeR finished at %s\n"%(ctime()))
	sys.stdout.flush()

	return
	
# End of main

if __name__ == '__main__':
	main()

