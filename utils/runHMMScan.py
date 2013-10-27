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
python runHMMScan.py --help
"""

USAGE = \
"""Usage: %prog <required_parameters> [options]

for more information, run:
python runHMMScan.py --help

runHMMScan.py is a utility script that helps prepare the hmmscan files used in BinGeR
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
from optparse import OptionParser, OptionGroup
import multiprocessing as mp
from subprocess import PIPE, Popen, call
from Bio import SeqIO
import re
from time import time

def runProdigal(x):
	infile, prodigal = x
	cmd = [prodigal, "-a", infile+'.faa', "-i", infile, "-p", "meta"]
	try:
		retcode = call(cmd, stdout = open(infile+'.gbk','a'), stderr = open(infile+'.log', 'a'))
	except OSError:
		sys.stderr.write( "OSError: fatal error running prodigal.\n" )
		exit(0)
	except ValueError:
		sys.stderr.write( "ValueError: fatal error running prodigal.\n" )
		exit(0)
	os.remove(infile+'.gbk')
	os.remove(infile+'.log')
	
# End of runProdigal

def runHMMScan(x):
	infile, outfile, db, prodigal, hmmscan = x
	runProdigal([infile, prodigal])
	faa = open(infile + '.prodigal.faa', 'w')
	tempfaa = open(infile + '.faa', 'r')
	for record in SeqIO.parse(tempfaa, 'fasta'):
		tag = record.id
		seq = record.seq
		tag = re.sub('\s\#\s(\d+)\s\#\s(\d+)\s\#\s(\-?\d+)\s\#\s', '|\1-\2|\3|', tag)
		faa.write('>%s\n%s\n'%(tag, seq))
	faa.close()
	tempfaa.close()
	
	os.remove(tempfaa)
	
	# run hmmscan here
	cmd = [hmmscan, "--tblout", outfile, "--cut_tc", db, faa]
	try:
		retcode = call(cmd, stdout = open(outfile+'.stdout', 'a'))
	except OSError:
		sys.stderr.write( "OSError: fatal error running hmmscan.\n" )
		exit(0)
	except ValueError:
		sys.stderr.write( "ValueError: fatal error running hmmscan.\n" )
		exit(0)
	os.remove(outfile+'.stdout')
	
def main(argv = sys.argv[1:]):

	parser = OptionParser(usage = USAGE, version="Version: " + __version__)
	
	# Required arguments
	
	requiredOptions = OptionGroup(parser, "Required options",
								"These options are required to run BinGeR, and may be supplied in any order.")
	
	requiredOptions.add_option("-i", "--infile", type = "string", metavar = "FILE",
							help = "Input fasta file.")

	requiredOptions.add_option("-o", "--outfile", type = "string", metavar = "DIR",
							help = "Output file with Z-scores of each input sequence.")

	requiredOptions.add_option("-d", "--hmm_db", type = "string", metavar = "FILE",
							help = "The single copy gene HMM model file.")

	parser.add_option_group(requiredOptions)

	# Optional arguments that need to be supplied if not the same as default
	optOptions = OptionGroup(parser, "Optional parameters",
						"There options are optional, and may be supplied in any order.")

	optOptions.add_option("-t", "--num_proc", type = "int", default = 1, metavar = 'INT',
							help = "Number of processors to use [default: 1].")
							
	optOptions.add_option("-p", "--prodigal", type = "string", default = "prodigal", metavar = 'STRING',
							help = "The prodigal program (Hyatt D. et al, Bioinformatics, 2010).")
							
	optOptions.add_option("-m", "--hmmscan", type = "string", default = "hmmscan", metavar = 'STRING',
							help = "The hmmscan program (http://hmmer.janelia.org).")
						
	parser.add_option_group(optOptions)
	
	(options, args) = parser.parse_args(argv)
	
	if options.infile is None:
		parser.error("An infile in fastA format is required!")
		exit(0)
	
	if not os.path.exists(options.infile):
		parser.error("Cannot find the infile you supplied: %s" % options.infile)
		exit(0)
		
	if options.num_proc < 1:
		parser.error("Number of processors must be integer >= 1, you supplied %i" % options.num_proc)
		exit(0)
		
	# test if prodigal exists
	prodigalTest = Popen(options.prodigal, shell=True, stderr=PIPE).stderr.read()
	if prodigalTest == None or len(prodigalTest) == 0:
		sys.stderr.write("FATAL: prodigal not found in path!\n")
		exit(0)
	
	# test if hmmscan exists
	hmmscanTest = Popen(options.hmmscan, shell=True, stdout=PIPE).stdout.read()
	if hmmscanTest == None or len(hmmscanTest) == 0:
		sys.stderr.write("FATAL: hmmscan not found in path!\n")
		exit(0)
		
	# run the main function here
	
	if options.num_proc == 1:
		runHMMScan([options.infile, options.outfile, options.hmm_db, options.prodigal, options.hmmscan])
	else:
		# split the input files
		tfhs = []
		infiles = []
		outfiles = []
		temp = str(time())
		for i in range(options.num_proc):
			tempfile = temp + '.fa.' + str(i+1)
			infiles.append(tempfile)
			outfiles.append(temp+'.hmmscan.'+str(i+1))
			tfh = open(tempfile, 'w')
			tfhs.append(tfh)
		
		ifh = open(options.infile, 'r')
		
		i = 0
		for record in SeqIO.parse(ifh, 'fasta'):
			i += 1
			index = i % options.num_proc
			tfhs[index].write('>%s\n%s\n' % (record.id, record.seq))
			
		ifh.close()
		for tfh in tfhs:
			tfh.close()
			
		# run mp jobs
		CMDs = [[infile, outfile, options.hmm_db, options.prodigal, options.hmmscan] \
					 for infile, outfile in zip(infiles, outfiles)]
		pool = mp.Pool(options.num_proc)
		pool.map_async(runHMMScan, CMDs)
		pool.close()
		pool.join()
		
		# merge results
		ofh = open(options.outfile, 'w')
		for outfile in outfiles:
			for line in open(outfile,'r'):
				ofh.write(line)
		ofh.close()
		for infile in files:
			os.remove(infile)
		for outfile in outfiles:
			os.remove(outfile)
	return
	
# End of main
	
if __name__ == '__main__':
	main()
	