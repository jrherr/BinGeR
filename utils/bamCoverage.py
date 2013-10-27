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
python bamCoverage.py --help
"""

USAGE = \
"""Usage: %prog <required_parameters> [options]
bamCoverage.py is a utility script that helps prepare the coverage files used in BinGeR
It is part of the BinGeR package under GNU 3 license. You are free to use and 
re-distribute it with the condition of keeping the script intact.

Additional information can also be found at:
https://github.com/luo-chengwei/BinGeR/wiki

If you use BinGeR in your work, please cite it as:
<BinGeR citation here>

Copyright: Chengwei Luo, Konstantinidis Lab, Georgia Institute of Technology, 2013
"""

import sys
import os
import glob
from optparse import OptionParser, OptionGroup
import multiprocessing as mp
import pysam

def coverageCal(x):
	bamFile, outfile = x
	samfh = pysam.Samfile(bamFile, 'rb')
	ofh = open(outfile, 'w')
	ofh.write('# Based on %s\n#ContigID\tContigLength\tNum_of_mapping_reads\tCoverage\n' % os.paht.abspath(bamFile))
	
	for contigID, contigLength in zip(samfh.references, samfh.lengths):
		numberReads = 0
		totalBPs = 0
		for read in samfh.fetch(contigID):
			if read.is_secondary:
				continue
			numberReads += 1
			totalBPs += read.rlen
		coverage = float(totalBPs)/contigLength
		ofh.write('%s\t%d\t%d\t%.4f\n'%(contigID, contigLength, numberReads, coverage))
		
	samfh.close()
	ofh.close()

# End of coverageCal

def main(argv = sys.argv[1:]):

	parser = OptionParser(usage = USAGE, version="Version: " + __version__)
	
	# Required arguments
	
	requiredOptions = OptionGroup(parser, "Required options",
								"These options are required to run BinGeR, and may be supplied in any order.")
	
	requiredOptions.add_option("-l", "--sample_list", type = "string", metavar = "FILE",
							help = "Text file containing all sample names, one per line")

	requiredOptions.add_option("-b", "--bams_dir", type = "string", metavar = "DIR",
							help = "Directory where all BAM files are")

	parser.add_option_group(requiredOptions)

	# Optional arguments that need to be supplied if not the same as default
	optOptions = OptionGroup(parser, "Optional parameters",
						"There options are optional, and may be supplied in any order.")

	optOptions.add_option("-c", "--coverage_dir", type = "string", default = "Coverage", metavar = "DIR",
							help = "Directory where coverage files are, naming follows \"sampleA.vs.sampleB.*.coverage\" convention. [Default: ./Coverage]")

	optOptions.add_option("-t", "--num_proc", type = "int", default = 1, metavar = 'INT',
							help = "Number of processor for BinGeR to use [default: 1].")
						

	parser.add_option_group(optOptions)
	
	(options, args) = parser.parse_args(argv)
	
	if options.sample_list is None:
		parser.error("A list of samples and a working directory are required!")
		exit(0)
		
	if options.bams_dir is None:
		parser.error("An BAM file directory is required to supply!")
		exit(0)
	
	if options.num_proc < 1:
		parser.error("Number of processors must be integer >= 1, you supplied %i" % options.num_proc)
		exit(0)
	
	samples = []
	try:
		sfh = open(options.sample_list, 'r')
		while 1:
			sample = sfh.readline().rstrip('\n')
			if not sample:
				break
			samples.append(sample)
		sfh.close()
	except:
		sys.stderr.write('FATAL: error in getting the list of all samples.\n')
		exit(0)
		
	# check if every file is there
	bamFiles = []
	covFiles = []
	for sampleA in samples:
		for sampleB in samples:
			files = glob.glob(options.bams_dir + '/' + sampleA + '.vs.' + sampleB + '.bam')
			if len(files) != 1:
				sys.stderr.write('FATAL: error in locating bam file for : %s vs % s' \
							% (sampleA, sampleB))
				exit(0)
			else:
				bamFiles += files
			outfile = options.coverage_dir + '/' + sampleA + '.vs.' + sampleB + '.coverage'
			covFiles.append(outfile)
		
	if not os.path.exists(options.coverage_dir):
		try:
			os.mkdir(options.coverage_dir)
		except IOError:
			sys.stderr.write('FATAL: cannot create output directory: %s\n'% options.coverage_dir)
			exit(0)
	
		
	CMDs = [[bamfile, covfile] for bamfile, covfile in zip(bamFiles, covFiles)]
	sys.stdout.write('Starting calculating the coverage...\n')
	pool = mp.Pool(options.num_proc)
	pool.map_async(coverageCal, CMDs)
	pool.close()
	pool.join()
	sys.stdout.write('Done. Results are saved at \n%s\n' % options.coverage_dir)

#End of main
	
if __name__ == '__main__':
	main()