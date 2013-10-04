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
import re
from string import maketrans
from operator import itemgetter
from Bio import SeqIO, Seq
from subprocess import call
import cPickle
import networkx as nx
import taxonomy

def nodePhylo(index, G, tTree, projInfo, options):
	nuc_dir = projInfo.out_dir + '/genes.nuc.clustered/'
	prot_dir = projInfo.out_dir + '/genes.prot.clustered/'
	trans = maketrans('ACGTacgt', 'TGCAtgca')
	
	# create lookup table for hmmscan contigs
	hmm_dict = {}
	for geneName in projInfo.HMM:
		for tag in projInfo.HMM[geneName]:
			hmm_dict[tag] = geneName

	nodes = {}
	for node in G.nodes(data = True):
		if 'HMM' not in node[1]:
			continue
		contigID = node[0]
		sample = node[1]['sample']
		if sample not in nodes:
			nodes[sample] = []
		nodes[sample].append(node)
		
	# holds all involved sequences
	if not options.quiet:
		sys.stdout.write('[initCore %i] Now extracting core genes...\n'%(index+1))
	
	# write to temp fasta file
	tempfasta = projInfo.out_dir + '/initCores/initCore.' + str(index + 1) + '.fa'
	if os.path.exists(tempfasta):
		pass
	else:
		sequences = {}
		for sample in projInfo.samples:
			assembly = projInfo.getAssemblyFile(sample)
			if sample not in nodes:
				continue
			sequences[sample] = {}
			afh = open(assembly, 'r')
			contigDict = SeqIO.to_dict(SeqIO.parse(afh, "fasta"))
			for node in nodes[sample]:
				record = contigDict[node[0]]
				seq = record.seq
				sequences[sample][node[0]] = seq
			afh.close()	
	
		# extract genes.
		candidateSeqs = []
		for sample in nodes:
			for node in nodes[sample]:
				contigID = node[0]
				for gene in node[1]['HMM']:
					tag = gene[0]
					start, end = re.search('(\d+)\-(\d+)', gene[1]).group(1,2)
					start = int(start)
					end = int(end)
					strand = gene[2]
					geneName = hmm_dict[gene[0]]
					geneSeq = sequences[sample][contigID][start-1:end]
				
					if strand != '1':
						geneSeq = str(geneSeq[::-1]).translate(trans)
					else:
						geneSeq = str(geneSeq)
					
					newTag = contigID+'|'+gene[1]+'|'+strand+'|'+geneName
					protSeq = Seq.Seq(geneSeq).translate()
					
					candidateSeqs.append((newTag, protSeq))
		
		sys.stdout.write('[initCore %i] Writing fasta\n'%(index+1))
		writeTempFasta(tempfasta, candidateSeqs)
		sys.stdout.write('[initCore %i] Done.\n'%(index+1))
	
	# run blat
	blatfile = tempfasta.replace('fa', 'blat')
	logfile = tempfasta.replace('fa', 'log')
	phyloGraphFile = tempfasta.replace('fa', 'graph.cpickle')
	database = projInfo.out_dir + '/genes.prot.clustered/allGenes.fasta'
	
	if os.path.exists(blatfile):
		pass
	else:
		if not options.quiet:
			sys.stdout.write('[initCore %i] Running blat.\n'%(index+1))
		runBLAT(tempfasta, database, blatfile, logfile, options.blat)
		if not options.quiet:
			sys.stdout.write('[initCore %i] Done.\n'%(index+1))
	
	# dict with contig->taxonomy mapping
	# interpret the results
	if os.path.exists(phyloGraphFile):
		if not options.quiet:
			sys.stdout.write('[initCore %i] Unpickling phyloGraph...\n'%(index+1))
		try:
			phfh = open(phyloGraph, 'rb')
			phyloGraph = cPickle.load(phfh)
			phfh.close()
		except:
			sys.stderr.write('FATAL: failure in unserializing phyloGraph.\n')
			exit(0)
		if not options.quiet:
			sys.stdout.write('[initCore %i] Done.\n'%(index+1))
	
	else:
		if not options.quiet:
			sys.stdout.write('[initCore %i] Rendering result...\n'%(index+1))
		
		phyloGraph = structPhylo(blatfile, tTree, projInfo)
		
		try:
			phfh = open(phyloGraphFile, 'wb')
			cPickle.dump(phyloGraph, phfh)
			phfh.close()
		except:
			sys.stderr.write('FATAL: failure in unserializing phyloGraph.\n')
			exit(0)
			
		if not options.quiet:
			sys.stdout.write('[initCore %i] Done.\n'%(index+1))
	
	# cleanup
#	os.remove(blatfile)
#	os.remove(logfile)
#	os.remove(tempfasta)

	return phyloGraph
	
# End of nodePhylo

def writeTempFasta(outfile, seqs):
	try:
		ofh = open(outfile, 'wb')
	except OSError:
		sys.stderr.write('FATAL: Cannot create temporary fasta file.\n')
	
	for seq in seqs:
		ofh.write('>'+seq[0]+'\n'+str(seq[1])+'\n')
	
	ofh.close()
	
# End of writeTempFasta

def runBLAT(query, db, outfile, logfile, blat):
	cmd = [blat, db, query, "-prot", "-out=blast8", "-noHead", outfile, logfile]
	try:
		retcode = call(cmd[:-1], stdout = open(cmd[-1],'a'))
	except OSError:
		sys.stderr.write( "OSError: fatal error running blat.\n" )
		exit(0)
	except ValueError:
		sys.stderr.write( "ValueError: fatal error running blat.\n" )
		exit(0)
# End of runBLAT

def structPhylo(blatfile, tTree, projInfo):
	## churn the blat file and pick the top 5 hits
	bfh = open(blatfile, 'r')
	blatRes = {}
	while 1:
		line = bfh.readline().rstrip('\n')
		if not line:
			break
		cols = line.split('\t')
		contig, start, end, strand, gene = re.search('(.+)\|(\d+)\-(\d+)\|(.+)\|(.+)$', cols[0]).group(1, 2, 3, 4, 5)
		taxid = int(re.search('(\d+)\.', cols[1]).group(1))
		start = int(start)
		end = int(end)
		geneLength = float(end - start + 1)/3
		identity = float(cols[2])
		percentageCov = float(cols[3])/geneLength
		bitscore = float(cols[-1])
		if bitscore < 50:
			continue
		
		if contig not in blatRes:
			blatRes[contig] = {}
		if gene not in blatRes[contig]:
			blatRes[contig][gene] = []
			
		if len(blatRes[contig][gene]) < 3:
			blatRes[contig][gene].append((taxid, percentageCov, identity, bitscore))
		elif len(blatRes[contig][gene]) == 3:
			bitscores = map(itemgetter(-1), blatRes[contig][gene])
			minScore = min(bitscores)
			if bitscore > minScore:
				# remove the lowest scored one.
				minIndex = bitscores.index(minScore)
				blatRes[contig][gene].pop(minIndex)
				blatRes[contig][gene].append((taxid, percentageCov, identity, bitscore))
	bfh.close()
	
	for contig in blatRes:
		for gene in blatRes[contig]:
			temp = [] 
			for hit in blatRes[contig][gene]:
				taxid = hit[0]
				taxonIDs, Ranks, sciNames = tTree.getTaxonomyPath(taxid)
				temp.append(taxonIDs)
			for i, taxonIDs in enumerate(temp):
				blatRes[contig][gene].insert(i*2, taxonIDs)
	
	# load partial tree/graph
	pTree = nx.Graph()
	edges = {}
	for contigID in blatRes:
		for gene in blatRes[contigID]:
			for i in range(0, len(blatRes[contigID][gene]), 2):
				taxonID = blatRes[contigID][gene][i]
				detail = blatRes[contigID][gene][i+1]
				bitscore = detail[-1]
				for i in range(len(taxonID)-1):
					nodeA = taxonID[i]
					if nodeA not in edges:
						edges[nodeA] = {}
					nodeB = taxonID[i+1]
					if nodeB not in edges[nodeA]:
						edges[nodeA][nodeB] = {'weight':0, 'genes':[]}
					edges[nodeA][nodeB]['weight'] += bitscore
					if i == 0:
						edges[nodeA][nodeB]['genes'].append((contigID, gene))
	e = []
	for nodeA in edges:
		for nodeB in edges[nodeA]:
			bitscore = edges[nodeA][nodeB]['weight']
			genes = edges[nodeA][nodeB]['genes']
			e.append((nodeA, nodeB, {'weight':bitscore, 'genes':genes}))
	pTree.add_edges_from(e)
	
	return pTree
	
# End of structPhylo
	
def weightedLCA(pTree, tTree):
	phylo = {}
	return phylo
	