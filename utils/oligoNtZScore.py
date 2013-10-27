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
python oligoNtZScore.py --help
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
from optparse import OptionParser, OptionGroup
import multiprocessing as mp
from Bio import SeqIO
import numpy as np
import scipy.stats
import string

def reverseComp(seq):
	seq=seq[-1::-1]
	seq=seq.translate(string.maketrans("ATGC","TACG"))
	return seq
# End of reverseComp

def tetraIndex(mer):
	base={}
	base['A']=0
	base['N']=0
	base['C']=1
	base['G']=2
	base['T']=3
	index=64*base[mer[0]]+16*base[mer[1]]+4*base[mer[2]]+base[mer[3]]
	return index
# End of tetraIndex

def triIndex(mer):
	base={}
	base['A']=0
	base['N']=0
	base['C']=1
	base['G']=2
	base['T']=3
	index=16*base[mer[0]]+4*base[mer[1]]+base[mer[2]]
	return index
# End of triIndex

def tetra(seq):
	freq=np.zeros(256)
	for i in range(len(seq)-3):
		tetraMer=seq[i:i+4]
		reTetraMer=reverseComp(tetraMer)
		indexA=tetraIndex(tetraMer)
		indexB=tetraIndex(reTetraMer)
		freq[indexA]+=1
		freq[indexB]+=1
	freq=freq/sum(freq)
	return freq
# End of tetra

def tri(seq):
	freq=np.zeros(64)
	for i in range(len(seq)-2):
		triMer=seq[i:i+3]
		reTriMer=reverseComp(triMer)
		indexA=triIndex(triMer)
		indexB=triIndex(reTriMer)
		freq[indexA]+=1
		freq[indexB]+=1
	freq=freq/sum(freq)
	return freq
# End of tri
	
def oligoStat(x):
	infile, outfile = x
	ifh = open(infile, 'r')
	ofh = open(outfile, 'w')
	for record in SeqIO.parser(ifh, 'fasta'):
		tag = record.id
		seq = record.seq
		tetraFreq = tetra(seq)
		triFreq = tri(seq)
		triZ = scipy.stats.zscore(triFreq)
		tetraZ = scipy.stats.zscore(tetraFreq)
		
		ofh.write('%s\t%d\n'%(tag, len(seq)))
		
		out=[]
		for f in triZ:
			out.append(str(f))
		triZLine="\t".join(out)
	
		out=[]
		for f in tetraZ:
			out.append(str(f))
		tetraZLine="\t".join(out)
	
		ofh.write(tag+'\n'+triZLine+'\n'+tetraZLine+'\n')
	ifh.close()
	ofh.close()
	
# End of oligoStat

def main(argv = sys.argv[1:]):

	parser = OptionParser(usage = USAGE, version="Version: " + __version__)
	
	# Required arguments
	
	requiredOptions = OptionGroup(parser, "Required options",
								"These options are required to run BinGeR, and may be supplied in any order.")
	
	requiredOptions.add_option("-i", "--infile", type = "string", metavar = "FILE",
							help = "Input fasta file.")

	requiredOptions.add_option("-o", "--outfile", type = "string", metavar = "DIR",
							help = "Output file with Z-scores of each input sequence.")

	parser.add_option_group(requiredOptions)

	# Optional arguments that need to be supplied if not the same as default
	optOptions = OptionGroup(parser, "Optional parameters",
						"There options are optional, and may be supplied in any order.")

	optOptions.add_option("-t", "--num_proc", type = "int", default = 1, metavar = 'INT',
							help = "Number of processors to use [default: 1].")
						
	parser.add_option_group(optOptions)
	
	(options, args) = parser.parse_args(argv)
	
	if not os.path.exists(options.infile):
		parser.error("Cannot find the infile you supplied: %s" % options.infile)
		exit(0)
		
	if options.num_proc < 1:
		parser.error("Number of processors must be integer >= 1, you supplied %i" % options.num_proc)
		exit(0)
		
	if options.num_proc == 1:
		oligoStat([options.infile, options.outfile])
	else:
		# split the input files
		tfhs = []
		infiles = []
		outfiles = []
		for i in range(options.num_proc):
			tempfile = 'temp.fa.' + str(i+1)
			infiles.append(tempfile)
			outfiles.append('temp.zscore.'+str(i+1))
			tfh = open(tempfile, 'w')
			tfhs.append(tfh)
		
		ifh = open(options.infile, 'w')
		
		i = 0
		for record in SeqIO.parse(ifh, 'fasta'):
			i += 1
			index = i % options.num_proc
			tfhs.write('>%s\n%s\n' % (record.id, record.seq))
			
		ifh.close()
		for tfh in tfhs:
			tfh.close()
			
		# run mp jobs
		CMDs = [[infile, outfile] for infile, outfile in zip(infiles, outfiles)]
		pool = mp.Pool(options.num_proc)
		pool.map_async(oligoStat, CMDs)
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
	
#End of main
	
if __name__ == '__main__':
	main()
	