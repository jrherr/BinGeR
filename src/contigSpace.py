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
import glob
import tarfile
import shutil
from subprocess import Popen, PIPE, call
import multiprocessing as mp
import cPickle
import networkx as nx
from Bio import SeqIO
import pysam
import numpy as np
from scipy.spatial import distance
from scipy.stats import norm

import phylo
from taxonomy import TaxonTree
from commPageRank import commPageRank


######################## CLASSES ######################
class ContigSpace(nx.Graph):
	def __init__(self):
		self.graphs = {}
		self.graph = nx.Graph()
		self.cores = {}
		self.bridges = nx.Graph()
		self.coreLookup = {}
		self.nodeLookup = {}
		self.candidateComp = []
	
	def __init__(self, samples):
		self.graphs = {}
		self.graph = nx.Graph()
		for sample in samples:
			self.graphs[sample] = nx.Graph()
		self.bridges = nx.Graph()
		self.coreLookup = {}
		self.nodeLookup = {}
		self.cores = {}
		self.candidateComp = []
	
	
	######################## UTILITY FUNCTIONS ##########################
	
	def getClusterName(self, x):
		try:
			return self.coreLookup[x][0]
		except KeyError, IndexError:
			sys.stderr.write('Error: cannot find the clusterID for %x'%x)
			exit(0)
	
	# End of getClusterName
			
	def getContigLength(self, x):
		try:
			return self.coreLookup[x][1]
			sys.stderr.write('Error: cannot find the length for %s\n'%x)
			exit(0)
		except KeyError, IndexError:
			sys.stderr.write('Error: cannot find the length for %s\n'%x)
	
	# End of getContigLength
	
	def getClusterLength(self, clusterID):
		L = 0
		for node in self.nodeLookup[clusterID]:
			L += node[1]
		return L
	
	# End of getClusterLength
	
	def maxRepLength(self, nodes):
		RepLength = {}
		for node in nodes:
			sample, clusterIDstr = re.search('(.+)\.(\d+)$', node).group(1, 2)
			clusterLength = self.getClusterLength(node)
			if sample not in RepLength:
				RepLength[sample] = 0
			RepLength[sample] += clusterLength
			
		return max(RepLength.values())
		
	# End of maxRepLength
	
	
	def addHMMScanInfoToSubgraph(self, sample, projInfo, options):
		
		# add hmmscan information
		if not options.quiet:
			sys.stdout.write('[%s] Adding HMMScan information...\n'%sample)
		hmmscanFile = projInfo.getHMMScanFile(sample)
		# read the hmmscan file, and load a lookup dict into hmmscanReg
		hmmscanReg = readHMMScanTable(hmmscanFile)
		nx.set_node_attributes(self.graphs[sample], 'HMM', hmmscanReg)
	
	# End of addHMMScanInfoToGraph
	
	def addEdgesFromCoverageFile(self, sample, projInfo, options):

		# add coverage information
		if not options.quiet:
			sys.stdout.write('[%s] Adding edges from coverage profile...\n' % sample)
	
		contigs = []
		data = []
	
		# go through all the coverage files that are against this assembly
		for index, readSample in enumerate(projInfo.samples):
			coverageFile = projInfo.getCoverageFile(readSample, sample)
			cfh = open(coverageFile, 'r')
	
			# discard the header lines (2)
			for k in range(2):
				disLine = cfh.readline()
	
			data.append([])
			while 1:
				line=cfh.readline().rstrip('\n')
				if not line:
					break
				ele = line.split('\t')
				length = int(ele[1])
				if length < options.minCovLength:
					continue
				if index == 0:
					contigs.append(ele[0])
				data[index].append(float(ele[-1]))
			cfh.close()
		
		matrix = np.array(data).transpose()
	
		thr = options.minCovCorrceof

		edges = edgesFromCoverageClustering(matrix, contigs, thr)
		
		self.graphs[sample].add_edges_from(edges, weight = 3)
		
		if not options.quiet:
			sys.stdout.write('\t%i edges added.\n' % len(edges))
			
		edges = []
		
	# End of addEdgesFromCoverageFiles
	
	def addEdgesFromZScoreFile(self, sample, projInfo, options):
				
		# add ZScore information
		
		if not options.quiet:
			sys.stdout.write('[%s] Adding edges from ZScores...\n'%sample)
	
		contigs = []
		tri = []
		tetra = []
		
		zscoreFile = projInfo.getZScoreFile(sample)
		zfh = open(zscoreFile, 'r')
		
		while 1:
			line=zfh.readline().rstrip('\n')
			if not line:
				break
			ele = line.split('\t')
			length = int(ele[1])
			contig = ele[0]
			
			if length < options.minZLength:
				disLine = zfh.readline()
				disLine = zfh.readline()
				continue
				
			contigs.append(ele[0])
			
			triLine = zfh.readline().rstrip('\n')
			triEle = triLine.split('\t')
			for i in range(len(triEle)):
				triEle[i] = float(triEle[i])
			tri.append(triEle)
			
			tetraLine = zfh.readline().rstrip('\n')
			tetraEle = tetraLine.split('\t')
			for i in range(len(tetraEle)):
				tetraEle[i] = float(tetraEle[i])
			tetra.append(tetraEle)
		zfh.close()
			
		thr = options.minZScoreCorrceof
	
		edges = edgesFromZScoreClustering(tri, tetra, contigs, thr)
		
		for e in edges:
			if self.graphs[sample].has_edge(e[0], e[1]):
				self.graphs[sample].add_edge(e[0], e[1], weight = 3.5)
			else:
				self.graphs[sample].add_edge(e[0], e[1], weight = 1.5)
		
		if not options.quiet:
			sys.stdout.write('\t%i edges added.\n' % len(edges))
			
		edges = []
		
	# End of addEdgesFromZScoreFile

	def loadGraphFromPickle(self, pickleFile, options):
		try:
			if not options.quiet:
				sys.stdout.write('Unserializing Graph...\n')
			cpfh = open(pickleFile, 'rb')
			self.graph = cPickle.load(cpfh)
			cpfh.close()
		except OSError:
			sys.stderr.write('FATAL: OSError in initializing graph from serialized subgraph!\n')
			sys.stderr.write('\tPlease check if you have the follow gpickle:\n\t%s\n'%pickleFile)
			exit(0)
		except cPickle.UnpicklingError:
			sys.stderr.write('FATAL: UnpicklingError in initializing graph from serialized subgraph!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: error in initializing contig space from serialized graph!\n')
			exit(0)
		
		if not options.quiet:
			sys.stdout.write('Done.\n')
	
	# End of loadGraphFromPickle
	
	def loadCoresFromPickle(self, pickleFile, options):
		try:
			if not options.quiet:
				sys.stdout.write('Unserializing cores...\n')
			cpfh = open(pickleFile, 'rb')
			self.cores = cPickle.load(cpfh)
			cpfh.close()
		except OSError:
			sys.stderr.write('FATAL: OSError in initializing cores from serialized subgraph!\n')
			sys.stderr.write('\tPlease check if you have the follow gpickle:\n\t%s\n'%pickleFile)
			exit(0)
		except cPickle.UnpicklingError:
			sys.stderr.write('FATAL: UnpicklingError in initializing cores from serialized core lookup!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: error in initializing cores from serialized core lookup!\n')
			exit(0)
		
		if not options.quiet:
			sys.stdout.write('Done.\n')
	
	# End of loadCoresFromPickle
	
	def loadCoreLookupFromPickle(self, pickleFile, options):
		try:
			if not options.quiet:
				sys.stdout.write('Unserializing core lookup...\n')
			cpfh = open(pickleFile, 'rb')
			self.coreLookup = cPickle.load(cpfh)
			cpfh.close()
		except OSError:
			sys.stderr.write('FATAL: OSError in initializing core lookup from serialized subgraph!\n')
			sys.stderr.write('\tPlease check if you have the follow gpickle:\n\t%s\n'%pickleFile)
			exit(0)
		except cPickle.UnpicklingError:
			sys.stderr.write('FATAL: UnpicklingError in initializing graph from serialized core lookup!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: error in initializing contig space from serialized core lookup!\n')
			exit(0)
		
		if not options.quiet:
			sys.stdout.write('Done.\n')
	
	# End of loadCoreLookupFromPickle
	
	def loadNodeLookupFromPickle(self, pickleFile, options):
		try:
			if not options.quiet:
				sys.stdout.write('Unserializing node lookup...\n')
			cpfh = open(pickleFile, 'rb')
			self.nodeLookup = cPickle.load(cpfh)
			cpfh.close()
		except OSError:
			sys.stderr.write('FATAL: OSError in initializing node lookup from serialized subgraph!\n')
			sys.stderr.write('\tPlease check if you have the follow gpickle:\n\t%s\n'%pickleFile)
			exit(0)
		except cPickle.UnpicklingError:
			sys.stderr.write('FATAL: UnpicklingError in initializing graph from serialized node lookup!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: error in initializing contig space from serialized node lookup!\n')
			exit(0)
		
		if not options.quiet:
			sys.stdout.write('Done.\n')
	
	# End of loadCoreLookupFromPickle
	
	def loadSubgraphsFromPickle(self, subgraphPickle):
		subgraphs = {}
		
		try:
			cpfh = open(subgraphPickle, 'rb')
			subgraphs = cPickle.load(cpfh)
			cpfh.close()
		except OSError:
			sys.stderr.write('FATAL: OSError in initializing node lookup from serialized subgraph!\n')
			sys.stderr.write('\tPlease check if you have the follow gpickle:\n\t%s\n'%subgraphPickle)
			exit(0)
		except cPickle.UnpicklingError:
			sys.stderr.write('FATAL: UnpicklingError in initializing graph from serialized node lookup!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: error in initializing contig space from serialized node lookup!\n')
			exit(0)
			
		return subgraphs
		
	# End of loadSubgraphsFromPickle
	
	def pickleGraph(self, pickleFile, options):
		# save the contig space if needed
		try:
			if not options.quiet:
				sys.stdout.write('Now pickling graph...\n')
			cpfh = open(pickleFile, 'wb')
			cPickle.dump(self.graph, cpfh)
			cpfh.close()
			if not options.quiet:
				sys.stdout.write('Contig space pickled, stored at:\n\t %s\n' % pickleFile)
		except OSError:
			sys.stderr.write('FATAL: OSError in serializing graph!\n')
			sys.stderr.write('\tPlease check if you have the privilege to generate this file:\n\t%s\n' % pickleFile)
			exit(0)
		except cPickle.PicklingError:
			sys.stderr.write('FATAL: PicklingError in serializing graph!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: Error in serializing graph!\n')
			exit(0)
		
		if not options.quiet:
			sys.stdout.write('Done.\n')
	
	# End of pickleGraph
	
	def pickleCores(self, pickleFile, options):
		try:
			if not options.quiet:
				sys.stdout.write('Now pickling cores...\n')
			cpfh = open(pickleFile, 'wb')
			cPickle.dump(self.cores, cpfh)
			cpfh.close()
			if not options.quiet:
				sys.stdout.write('Cores pickled, stored at:\n\t %s\n' % pickleFile)
		except OSError:
			sys.stderr.write('FATAL: OSError in serializing cores!\n')
			sys.stderr.write('\tPlease check if you have the privilege to generate this file:\n\t%s\n' % pickleFile)
			exit(0)
		except cPickle.PicklingError:
			sys.stderr.write('FATAL: PicklingError in serializing cores!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: Error in serializing cores!\n')
			exit(0)
		
		if not options.quiet:
			sys.stdout.write('Done.\n')
	
	# End of pickleCores
	
	def pickleCoreLookup(self, pickleFile, options):
		# save the contig space if needed
		try:
			if not options.quiet:
				sys.stdout.write('Now pickling graph...\n')
			cpfh = open(pickleFile, 'wb')
			cPickle.dump(self.coreLookup, cpfh)
			cpfh.close()
			if not options.quiet:
				sys.stdout.write('Core lookup pickled, stored at:\n\t %s\n' % pickleFile)
		except OSError:
			sys.stderr.write('FATAL: OSError in serializing graph!\n')
			sys.stderr.write('\tPlease check if you have the privilege to generate this file:\n\t%s\n' % pickleFile)
			exit(0)
		except cPickle.PicklingError:
			sys.stderr.write('FATAL: PicklingError in serializing core lookup!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: Error in serializing core lookup!\n')
			exit(0)
		
		if not options.quiet:
			sys.stdout.write('Done.\n')
	
	# End of pickleCoreLookup
	
	def pickleNodeLookup(self, pickleFile, options):
		# save the contig space if needed
		try:
			if not options.quiet:
				sys.stdout.write('Now pickling graph...\n')
			cpfh = open(pickleFile, 'wb')
			cPickle.dump(self.nodeLookup, cpfh)
			cpfh.close()
			if not options.quiet:
				sys.stdout.write('Node lookup pickled, stored at:\n\t %s\n' % pickleFile)
		except OSError:
			sys.stderr.write('FATAL: OSError in serializing graph!\n')
			sys.stderr.write('\tPlease check if you have the privilege to generate this file:\n\t%s\n' % pickleFile)
			exit(0)
		except cPickle.PicklingError:
			sys.stderr.write('FATAL: PicklingError in serializing node lookup!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: Error in serializing node lookup!\n')
			exit(0)
		
		if not options.quiet:
			sys.stdout.write('Done.\n')
	
	
	############### Major functions that initiates graph #################
	def initSubgraphs(self, projInfo, options):
		
		if not options.quiet:
			sys.stdout.write('Constructing contig space for each sample...\n')
		
		subgraphDir = projInfo.out_dir + '/subgraphs/'
		
		if not os.path.exists(subgraphDir):
			try:
				os.mkdir(subgraphDir)
			except OSError:
				sys.stderr.write('FATAL [OSError]: Error in creating subgraphs dir, please check if you have privilege at %s.\n'%projInfo.out_dir)
				exit(0)
			except:
				sys.stderr.write('FATAL: Error in creating subgraphs dir, please check if you have privilege at %s.\n'%projInfo.out_dir)
				exit(0)
		
		coreLookupPickle = projInfo.out_dir + '/coreLookup.init.cpickle'
		nodeLookupPickle = projInfo.out_dir + '/nodeLookup.init.cpickle'
		graphPickle = projInfo.out_dir + '/graph.init.cpickle'
		
		if os.path.exists(coreLookupPickle) and os.path.exists(nodeLookupPickle) and os.path.exists(graphPickle):
			if not options.quiet:
				sys.stdout.write('Start loading information...\n')
			self.loadCoreLookupFromPickle(coreLookupPickle, options)
			self.loadNodeLookupFromPickle(nodeLookupPickle, options)
			self.loadGraphFromPickle(graphPickle, options)
			
		else:
			for sample in projInfo.samples:
				self.loadSubgraphs(sample, projInfo, options)
		
			self.pickleGraph(graphPickle, options)
			self.pickleCoreLookup(coreLookupPickle, options)
			self.pickleNodeLookup(nodeLookupPickle, options)
		
		sys.stdout.write('Done initializing the contig space.\n')
	
	# End of initSubgraphs
	
	def loadSubgraphs(self, sample, projInfo, options):
		# init the nodes first
		assemblyFile = projInfo.getAssemblyFile(sample)
		afh = open(assemblyFile, 'r')
		contigs = []
		for seq in SeqIO.parse(afh, "fasta"):
			tag = seq.id.split(' ')[0]
			length = len(seq.seq)
			contigs.append((tag,{'length':length, 'sample':sample}))
		
		afh.close()
		
		self.graphs[sample].add_nodes_from(contigs)
		sys.stdout.write('[%s] %i contigs added.\n'%(sample, len(contigs)))
			
		# add coverage linkages
		self.addEdgesFromCoverageFile(sample, projInfo, options)
		
		# add ZScore linkages
		self.addEdgesFromZScoreFile(sample, projInfo, options)
			
		# add HMMScan results
		self.addHMMScanInfoToSubgraph(sample, projInfo, options)
			
		# summary
		numNodes = nx.number_of_nodes(self.graphs[sample])
		numEdges = nx.number_of_edges(self.graphs[sample])
			
		if not options.quiet:
			sys.stdout.write('[%s] Subgraph has %i nodes and %i edges.\n'%(sample, numNodes, numEdges))
			
		# cluster and store the subgraphs
		self.shrinkSubgraph(sample, projInfo, options)
			
		# summary
		numNodes = nx.number_of_nodes(self.graph)
			
		if not options.quiet:
			sys.stdout.write('Graph has %i nodes after adding %s.\n\n'%(numNodes, sample))
		
	# End of loadSubgraphs
	
	def shrinkSubgraph(self, sample, projInfo, options):
		## Now collapse it, and make coreLookup and nodeLookup
		## coreLookup is a contigID->(clusterName, length) mapping
		## and nodeLookup is a clusterName->tuples((contigID, length)...)
		## of cluster (list) mapping
		
		subgraphs = {}
		clusterID = 0
		clusterNames = []
		
		if not options.quiet:
			sys.stdout.write('[%s] Now shrinking the subgraph and constructing lookup tables...\n'%sample)
		
		for subgraph in nx.connected_component_subgraphs(self.graphs[sample]):
			clusterID += 1
			clusterName = sample + '.' + str((clusterID))
			clusterNames.append(clusterName)
			self.nodeLookup[clusterName] = []
			subgraphs[clusterName] = subgraph.copy()
			# build up lookups
			nodes = subgraph.nodes(data = True)
			for node in nodes:
				contigID = node[0]
				sample = node[1]['sample']
				length = node[1]['length']
				self.nodeLookup[clusterName].append((contigID, length))
				self.coreLookup[contigID] = ((clusterName, length))
		
		self.graphs[sample].clear()
		self.graph.add_nodes_from(clusterNames)   # init self.graph with nodes(clusterNames)
			
		# save the contig space if needed
		pickleFile = projInfo.out_dir + '/subgraphs/' + sample + '.subgraphs.cpickle'
		try:
			if not options.quiet:
				sys.stdout.write('[%s] Now pickling the subgraphs...\n'%sample)
			cpfh = open(pickleFile, 'wb')
			cPickle.dump(subgraphs, cpfh)
			cpfh.close()
			if not options.quiet:
				sys.stdout.write('[%s] Subgraphs pickled.\n' % sample)
		except OSError:
			sys.stderr.write('FATAL: OSError in serializing subgraph!\n')
			sys.stderr.write('\tPlease check if you have the privilege to generate this file:\n\t%s\n' % pickleFile)
			exit(0)
		except cPickle.PicklingError:
			sys.stderr.write('FATAL: PicklingError in serializing subgraph!\n')
			exit(0)
		except:
			sys.stderr.write('FATAL: Error in serializing contig space!\n')
			exit(0)
		
	# End of shrinkSubgraph
	
	############### Major function to pick cores #############
	def forgeCores(self, projInfo, options):
		"""
		This function goes through the init graph formed in initGraph(),
		and it loades BLAT result to connect different samples, it then
		loads detailed information to pass to refineCores
		"""
		
		# if it is already pickled, then just load from there
		
		# bridgePicke -> self.bridge: nx.Graph(), nodes: (contigID, {'length':L})
		# edges: contig-contig edges, (contigID_A, contigID_B, {'weight':W, 'alignLength':L})
		
		bridgePickle = projInfo.out_dir + '/bridges.cpickle'
		
		# subgraphsPickle -> self.graphs
		# dict keyed by clusterName ('sample.clusterID')
		# values are nx.Graph(), nodes are (contigID, {'length':L, 'sample':sample})
		# for nodes with HMMScan results, attr dict is {'length':L, 'sample':sample, 'HMM':hmmscanReg}
		# where hmmscanReg is a dict: {'contigID':(hmmID, position, strand)}, 
		# postion: 'start-end', strand: '1' or '-1'
		# edges are (contigID_A, contigID_B, {'weight':W})
		
		subgraphsPickle = projInfo.out_dir + '/subgraphs.cpickle'
		
		# candidateComponentsPickle->self.candidateComp
		# dict keyed by sample, values are lists of int type clusterIDo
		
		candidateComponentsPickle = projInfo.out_dir + '/candidateComp.cpickle'
		
		if os.path.exists(bridgePickle) and os.path.exists(subgraphsPickle) and os.path.exists(candidateComponentsPickle):
			if not options.quiet:
				sys.stdout.write('Unserializing pickled data...\n')
			
			cpfh = open(bridgePickle, 'rb')
			try:
				self.bridges = cPickle.load(cpfh)
				cpfh.close()
			except:
				sys.stderr.write('FATAL: Error occurred during unserializing bridges.\n')
				exit(0)
				
			cpfh = open(subgraphsPickle, 'rb')
			try:
				self.graphs = cPickle.load(cpfh)
				cpfh.close()
			except:
				sys.stderr.write('FATAL: Error occurred during unserializing subgraphs.\n')
				exit(0)
				
			cpfh = open(candidateComponentsPickle, 'rb')
			try:
				self.candidateComp = cPickle.load(cpfh)
				cpfh.close()
			except:
				sys.stderr.write('FATAL: Error occurred during unserializing candidate components.\n')
				exit(0)
			
			if not options.quiet:
				sys.stdout.write('Done.\n')
				
		
		else:
			# load BLAT edges first
			if not options.quiet:
				sys.stdout.write('Loading blat edges...\n')
			
			blatEdgePickle = projInfo.out_dir + '/BlatEdges.cpickle'
			if not os.path.exists(blatEdgePickle):
				generateBlatEdgePickle(projInfo, options)
			blatEdges = loadBlatEdgesFromPickle(blatEdgePickle, options)
		
			# construct a temporary graph to clean up the blat edges
			if not options.quiet:
				sys.stdout.write('Done. Now cleaning the edges...\n')
		
		
			node_attributes = {}
			# add edges
			for sampleA in blatEdges:
				for sampleB in blatEdges[sampleA]:
					links = blatEdges[sampleA][sampleB]
					self.bridges.add_edges_from(links)
		
			# add nodes attribute(length) to the bridges
			for node in self.bridges.nodes():
				length = self.getContigLength(node)
				node_attributes[node] = length
			nx.set_node_attributes(self.bridges, 'length', node_attributes)
		
			"""
			iterate through all nodes, calculate the "thickness" of per-base
			pileup of every contig(node).
		
			Theoretically, the number shouldn't be larger than the number of samples - 1	
			if it is larger than # samples - 1, it is an indication that this contig will cause
			some chimerism.
		
			The promiscuity coefficient for the ith contig is defined as:
		
			P(i) = sigma(Ai)/Li, 
		
			where sigma(A) is the sum of the aligned length of all contigs that
			share an edge with the ith contig, and Li is the length in bp of the ith contig.
			"""
	
			if not options.quiet:
				sys.stdout.write('Now trimming promiscuous contigs...\n')
			
			for node in self.bridges.nodes(data = True):
				nodeID = node[0]
				nodeLength = node[1]['length']
				neighborLength = 0
				neighbors = self.bridges[nodeID]
				for e in neighbors:
					neighborLength += neighbors[e]['alignLength']
				promiscuity = float(neighborLength)/nodeLength		
				if promiscuity > 1.5 * len(projInfo.samples):
					self.bridges.remove_node(node)
		
			numSub = 0
			for c in nx.connected_components(self.bridges):
				if len(c) > 5:
					numSub+=1
			sys.stdout.write('%i connected components with 5+ members from blat graph.\n'%numSub)
		
			# Now from bridges, outputs the edges to be integrated, also
			# stores the edges to self.bridges
			blatEdges = []
			for edge in self.bridges.edges(data = True):
				clusterNameA = self.getClusterName(edge[0])
				clusterNameB = self.getClusterName(edge[1])
			
				blatEdges.append((clusterNameA, clusterNameB, {'weight':0.1}))
		
			# now add these blatEdges into self.graph
			self.graph.add_edges_from(blatEdges)
			sys.stdout.write('%i nodes and %i edges in self.graph\n'%(nx.number_of_nodes(self.graph), nx.number_of_edges(self.graph)))
		
			# now go through core candidates
			if not options.quiet:
				sys.stdout.write('\nNow forging init cores...\n')
		
			candidateClusters = {}        # keys would be the samples, values would be the clusterIDs in list
		
			for component in nx.connected_components(self.graph):
			
				maxRepLength = self.maxRepLength(component)
				if maxRepLength < 0.5 * options.min_core:
					continue
				
				# them as candidates for later
				# to expand the component,  cluster further with MCL, 
				# evaluate with HMMScan results, and form the core.
			
				self.candidateComp.append(component)

				for clusterName in component:
					sample, clusterIDstr = re.search('(.+)\.(\d+)$', clusterName).group(1, 2)
					clusterID = int(clusterIDstr)
					if sample not in candidateClusters:
						candidateClusters[sample] = []
					candidateClusters[sample].append(clusterID)
		
			if not options.quiet:
				for sample in projInfo.samples:
					try:
						sys.stdout.write('[%s] %i candidate clusters\n'%(sample, len(candidateClusters[sample])))
					except KeyError:
						sys.stdout.write('[%s] 0 candidate cluster\n'%sample)
					
			# retrieving subgraph information from pickled data
			self.graphs = {}     # holds all subgraphs keyed by pattern 'sample.clusterID'
			if not options.quiet:
				sys.stdout.write('\nLoading subgraph information...\n')
			
			for subgraphPickle in glob.glob(projInfo.out_dir + '/subgraphs/*.subgraphs.cpickle'):
				sample = re.search('.+\/(.+)\.subgraphs\.cpickle', subgraphPickle).group(1)
				if not options.quiet:
					sys.stdout.write('[%s] ' % sample)
				subgraphs = self.loadSubgraphsFromPickle(subgraphPickle)
			
				if not options.quiet:
					sys.stdout.write('Loaded %i subgraphs, '%len(subgraphs))

				for clusterName in subgraphs:
					clusterID = int(re.search('.+\.(\d+)', clusterName).group(1))
					if clusterID in candidateClusters[sample]:
						self.graphs[clusterName] = subgraphs[clusterName]
				if not options.quiet:
					sys.stdout.write('total number of passed subgraphs is now %i\n' % len(self.graphs))
				
			if not options.quiet:
				sys.stdout.write('Done.\n')
			
			
			# save to pickles (self.bridges, self.graphs, self.candidateComp)
			if not options.quiet:
				sys.stdout.write('\nNow serializing results...\n')
		
			cpfh = open(bridgePickle, 'wb')
			try:
				cPickle.dump(self.bridges, cpfh)
			except:
				sys.stderr.write('FATAL: Error occurred during pricking bridges.\n')
				exit(0)
			cpfh.close()
		
			cpfh = open(subgraphsPickle, 'wb')
			try:
				cPickle.dump(self.graphs, cpfh)
			except:
				sys.stderr.write('FATAL: Error occurred during pricking subgrpahs.\n')
				exit(0)
			cpfh.close()
		
			cpfh = open(candidateComponentsPickle, 'wb')
			try:
				cPickle.dump(self.candidateComp, cpfh)
			except:
				sys.stderr.write('FATAL: Error occurred during pricking candidate components.\n')
				exit(0)
			cpfh.close()
		
			if not options.quiet:
				sys.stdout.write('Done.\n')
		
	# End of forgeCores
	
	############### Major function to refine cores #############
	def refineCores(self, projInfo, options):
		
		"""
		# after all edges being transformed, it goes through every connected component,
		# if in none of the samples the connected component exceeds the size defined
		# by options.min_core, then it won't be considered further.
		#
		# For the cores that passed the filtering, each core will first be assess by their
		# single-copy gene distribution, as well as the shape of the cluster. If necessary,
		# it will go through a custom Markov Clustering algorithm (Van Dongen S., 2000;
		# see micans.org/mcl/ for theoretical details) to split them until converge.
		#
		# The cores will be output to self.core attribute for later analysis
		"""
		
		if not options.quiet:
			sys.stdout.write('Now refining cores...\n')
		
		# results are stored at out_dir/initCore/
		
		initCoresPath = projInfo.out_dir + '/initCores'
		if not os.path.exists(initCoresPath):
			try:
				os.mkdir(initCoresPath)
			except OSError:
				sys.stderr.write('FATAL: error when creating init core path.\n')
				exit(0)
				
		# if init cores exists, then from there
		initCoresPickles = glob.glob(initCoresPath + '/*.cpickle')
		initCores = []
		if len(initCoresPickles) > 0:
			if not options.quiet:
				sys.stdout.write('Loading init cores from pickle...\n')
			
			for index, initCoresPickle in enumerate(initCoresPickles):
				if index == 0:
					continue
				initCore = nx.Graph()
				try:
					cpfh = open(initCoresPickle, 'rb')
					initCore = cPickle.load(cpfh)
					cpfh.close()
					initCores.append(initCore)
				except:
					sys.stderr.write('FATAL: Error during unpickling init cores.\n')
					exit(0)
			if not options.quiet:
				sys.stdout.write('Done.\n')
				
		else:
			for index, component in enumerate(self.candidateComp):
				if not options.quiet:
					sys.stdout.write('[initCore %i] Graph initialization.\n'%(index+1))
				# expand the component.
				coreGraph = nx.Graph()
				for clusterName in component:
					subgraph = self.graphs[clusterName]
					coreGraph.add_nodes_from(subgraph.nodes(data = True))
					coreGraph.add_edges_from(subgraph.edges(data = True))
				
				# update BLAT edges from self.bridges
				blatEdges = []
				for edge in self.bridges.edges(data = True):
					nodeA = edge[0]
					nodeB = edge[1]
					if coreGraph.has_node(nodeA) and coreGraph.has_node(nodeB):
						weight = edge[2]['weight']
						blatEdges.append((nodeA, nodeB, {'weight':weight}))
				
				coreGraph.add_edges_from(blatEdges)
				self.bridges.remove_edges_from(blatEdges)
				
				if not options.quiet:
					sys.stdout.write('Pickling init core...\n')
				
				initCores.append(coreGraph)
				
				# pickle the initCores
				initCoresPickle = initCoresPath + '/initCore.' + str(index+1) + '.cpickle'
				try:
					cpfh = open(initCoresPickle, 'wb')
					cPickle.dump(coreGraph, cpfh)
					cpfh.close()
				except:
					sys.stderr.write('FATAL: Error during pickling init cores.\n')
					exit(0)
		
		# go through each initCore and apply community PageRank
		
		if not options.quiet:
			sys.stdout.write('Now loading taxonomy tree...\n')
		tTree = TaxonTree()
		tTree.loadTreeFromNodeLib(projInfo.DBs['ncbiNodes'], projInfo.DBs['ncbiSciNames'])
		if not options.quiet:
			sys.stdout.write('Done.\n')
		
		# extract db files to out_dir
		if not options.quiet:
			sys.stdout.write('Now preparing single copy genes...\n')
		try:
			nucTar = tarfile.open(projInfo.DBs['nuc'], 'r')
			protTar = tarfile.open(projInfo.DBs['prot'], 'r')
		except:
			sys.stderr.write('FATAL: failure in opening tarfile.\n')
			exit(0)
		
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
		
		# create merged files
		catNuc = open(projInfo.out_dir + '/genes.nuc.clustered/allGenes.fasta', 'w')
		catProt = open(projInfo.out_dir + '/genes.prot.clustered/allGenes.fasta', 'w')
		
		for file in glob.glob(projInfo.out_dir + '/genes.nuc.clustered/*.fa'):
			shutil.copyfileobj(open(file, 'rb'), catNuc)
		catNuc.close()
		
		for file in glob.glob(projInfo.out_dir + '/genes.prot.clustered/*.fa'):
			shutil.copyfileobj(open(file, 'rb'), catProt)
		catProt.close()
		
		if not options.quiet:
			sys.stdout.write('Done.\n')
			
		# go through each initCore and refine them.
		coreID = 0
		for coreIndex, initCore in enumerate(initCores):
			# get the taxonomy affiliation
			# returned seed nodes is a dict keyed by weighted LCA and valued by list of contigIDs
			
			pTree = phylo.nodePhylo(coreIndex, initCore, tTree, projInfo, options)
			
			continue
			
			seedNodes = phylo.weightedLCA(pTree)
			
			# refine the graph using community PageRank
			sys.stdout.write('Running personalized PageRank algorithms now.\n')
			refinedCoreGraph = commPageRank(initCore, seedNodes, options)
			
			# output to self.cores
			# need to output taxonomy affiliation as well.
			for core in nx.connected_components(refinedCoreGraph):
				coreID += 1
				self.cores[coreID] = []
				for contigID in core:
					self.cores[coreID].append(contigID)
	
				if not options.quiet:
					sys.stdout.write('initCore %i has %i nodes \n'%(coreID, len(core)))
		
		return
		
		# clean up self.graphs to release RAM
		for clusterName in self.graphs:
			self.graphs[clusterName].clear()
			
		# pickle self.cores
		self.pickleCores(corePickle, options)
		
		if not options.quiet:
			sys.stdout.write('Done.\n')

		# clean up the directory
		try:
			shutil.rmtree(projInfo.out_dir + '/genes.nuc.clustered/')
		except:
			sys.stderr.write('WARNING: Failure in removing intermediate nucleotide sequences directory.\n')
			
		try:
			shutil.rmtree(projInfo.out_dir + '/genes.prot.clustered/')
		except:
			sys.stderr.write('WARNING: Failure in removing intermediate protein sequences directory.\n')	
	# End of refineCores


	############ major function to recruit contigs to cores #############
	def recruitContigs(self, projInfo, options):
		print 'Code under construction\n'

	# End of recruitContigs

	def outputBins(self, projInfo, options):
		print 'Code under construction\n'

	# End of outputBins

	def extractReadsForBins(self, projInfo, options):
		print 'Code under construction\n'

	# End of extractReadsForBins
	
# End of class ContigSpace

######################## FUNCTIONS ####################

def runBLAT(cmd):
	try:
		retcode = call(cmd[:-1], stdout = open(cmd[-1],'a'))
	except OSError:
		sys.stderr.write( "OSError: fatal error running blat.\n" )
		exit(0)
	except ValueError:
		sys.stderr.write( "ValueError: fatal error running blat.\n" )
		exit(0)
        
# End of runBLAT

def readHMMScanTable(hmmscanFile):
	try:
		hmmfh = open(hmmscanFile, 'r')
	except IOError:
		sys.stderr.write('FATAL: cannot open HMMScan file: %s\n'%hmmscanFile)
		exit(0)
	
	hmmscanReg = {}
	while 1:
		line = hmmfh.readline().rstrip('\n')
		if not line:
			break
		if line[0] == '#':
			continue
		col = re.split('\s+', line)
		hmmID = col[0]
		contig = col[2].split('_')[0]
		eles = col[2].split('|')
		position = eles[1]
		strand = eles[2]
		
		if contig not in hmmscanReg:
			hmmscanReg[contig] = []
		hmmscanReg[contig].append((hmmID, position, strand))
	hmmfh.close()
	
	return hmmscanReg

# End of readHMMScanTable

def generateBlatEdgePickle(projInfo, options):
		
	blat = options.blat
	num_proc = options.num_proc
	save_intermediates = options.save_intermediates
	quiet = options.quiet

	###################### BLAT SEARCH #######################
	# run blat to find out how to group contigs together
	if not quiet:
		sys.stdout.write('Start conglomerate contigs using blat!\n')
	combinations = []
	for i, sampleA in enumerate(projInfo.samples):
		for j, sampleB in enumerate(projInfo.samples):
			if i >= j:
				continue
			combinations.append((sampleA, sampleB))
	
	# run blat here
	if not os.path.exists(projInfo.out_dir+'/blat'):
		try:
			os.mkdir(projInfo.out_dir+'/blat')
		except OSError:
			os.stderr.write("FATAL: cannot create blat in output directory: %s, please check your previlege\n"%projInfo.out_dir)
			exit(0)
	
	blatDir = os.path.abspath(projInfo.out_dir+'/blat')
	blatCMDs = []

	for combination in combinations:	
		outfile = blatDir+'/'+combination[0]+'.vs.'+combination[1]+'.blat'
		if os.path.exists(outfile):
			continue
		logfile = blatDir+'/'+combination[0]+'.vs.'+combination[1]+'.log'
		query = projInfo.getAssemblyFile(combination[0])
		db = projInfo.getAssemblyFile(combination[1])
		blatCMD = [blat, db, query, "-out=psl", "-noHead", "-fastMap", "-minIdentity=95", outfile, logfile]
		blatCMDs.append(blatCMD)

	if len(blatCMDs) > 0:	
		pool = mp.Pool(num_proc)
		rval = pool.map_async(runBLAT, blatCMDs)
		pool.close()
		pool.join()

	if not quiet:
		sys.stdout.write('Blat search finished, now calculating edges...\n')
	# blat search done
	################################################################
	
	# read the blat output, figure out all the edges and the weight
	links = {}
	for blatFile in glob.glob(blatDir+'/*.blat'):
		sampleA, sampleB = re.search('.+\/(.+)\.vs\.(.+)\.blat', blatFile).group(1, 2)
		if sampleB not in links:
			links[sampleB] = {}
		if sampleA not in links[sampleB]:
			links[sampleB][sampleA] = []
			
		bfh = open(blatFile, 'r')
		while 1:
			line = bfh.readline().rstrip('\n')
			if not line:
				break
			col = line.split('\t')
			matches = int(col[0])
			mismatches = int(col[1])
			qBaseInsert = int(col[6])
			identity=float(matches)/(matches + mismatches + qBaseInsert)
			if identity < 0.98 or matches < 400:
				continue
			qName = col[9]
			tName = col[13]
			qSize = int(col[10])
			tSize = int(col[14])
			qStart, qEnd = (int(col[11]), int(col[12]))
			tStart, tEnd = (int(col[15]), int(col[16]))
		
			if qStart > qEnd :
				qStart, qEnd = qEnd, qStart
			if tStart > tEnd :
				tStart, tEnd = tEnd, tStart
		
			alignedPercentage = 100.0*(qEnd - qStart + 1)/min(qSize, tSize)
			
			if alignedPercentage >= 90:      # if aligned > 90% of the shorter contigs, pass
				alignLength = qEnd - qStart + 1
				weight = 1
				if alignLength >= 500 and alignLength < 1500:
					weight = 2
				elif alignLength >=1500 and alignLength < 3000:
					weight = 2.5
				elif alignLength >= 3000:
					weight = 3
				links[sampleB][sampleA].append((qName, tName, {'weight': weight, 'alignLength':alignLength}))
			else:                            # otherwise, it could be sideways alignment
				q5Edge = qStart - 1
				q3Edge = qSize - qEnd + 1
			
				t5Edge = tStart - 1
				t3Edge = tSize - tStart + 1
			
				if max(q3Edge, t5Edge) < 100 or max(q5Edge, t3Edge) < 100:
					alignLength = qEnd - qStart + 1
					weight = 1.5
					if alignLength >= 500 and alignLength < 1500:
						weight = 2
					elif alignLength >=1500 and alignLength < 3000:
						weight = 2.5
					elif alignLength >= 3000:
						weight = 3
					links[sampleB][sampleA].append((qName, tName, {'weight': weight, 'alignLength':alignLength}))	

		bfh.close()
			
	pickleFile = projInfo.out_dir + '/BlatEdges.cpickle'
	cpfh = open(pickleFile, 'wb')
	if not options.quiet:
		sys.stdout.write('Now pickling all the BLAT edges...\n')
	cPickle.dump(links, cpfh)
		
	if not options.quiet:
		sys.stdout.write('Done.\n')
	cpfh.close()
		
	# End of generateBlatEdgePickle

def loadBlatEdgesFromPickle(blatPickle, options):
	try:
		if not options.quiet:
			sys.stdout.write('Unserializing BLAT edges...\n')
		cpfh = open(blatPickle, 'rb')
		edges = cPickle.load(cpfh)
		cpfh.close()
	except OSError:
		sys.stderr.write('FATAL: OSError in initializing contig space from serialized subgraph!\n')
		sys.stderr.write('\tPlease check if you have the follow gpickle:\n\t%s\n'%pickleFile)
		exit(0)
	except cPickle.UnpicklingError:
		sys.stderr.write('FATAL: UnpicklingError in initializing contig space from serialized subgraph!\n')
		exit(0)
#	except:
#		sys.stderr.write('FATAL: error in initializing contig space from serialized subgraph!\n')
#		exit(0)
		
	if not options.quiet:
		sys.stdout.write('Done.\n')
		
	return edges
		
# End of loadBlatEdgesFromPickle

def edgesFromCoverageClustering(data, columnLabels, threshold):
	edges = []
	# calculate the correlation coeffecient matrix
	corrcoef = np.corrcoef(data)
	# select elements that passed corrcoef thr
	indexArray = np.where(corrcoef > threshold)
	sys.stdout.write('%i pre-edges...\n'%len(indexArray[0]))
	# calculate normed distance
	# form a new np.array first to reduce computational cost
	
	for index in range(len(indexArray[0])):
		i = indexArray[0][index]
		j = indexArray[1][index]
		normMinkowski = distance.euclidean(data[i,], data[j,])/np.sqrt(np.dot(data[i,], data[j,]))
		if normMinkowski < 0.1:
			edges.append((columnLabels[i], columnLabels[j]))
#	sys.stdout.write('corroef done!\n')
	
	return edges
	
# End of edgesFromCoverageClustering

def edgesFromZScoreClustering(tri, tetra, columnLabels, threshold):
	# calculate the correlation coeffecient matrices
	triCorrCoef = np.corrcoef(tri)
#	sys.stdout.write('tri-nt corroef done!\n')
	tetraCorrCoef = np.corrcoef(tetra)
#	sys.stdout.write('tetra-nt corroef done!\n')
	
	edges = []
	
	# make a boolean matrix
	indexArray = np.where((triCorrCoef > threshold) & (tetraCorrCoef > threshold))
	for i in range(len(indexArray[0])):
		x = indexArray[0][i]
		y = indexArray[1][i]
		edges.append((columnLabels[x], columnLabels[y]))
#	sys.stdout.write('edges addd\n')
	
	return edges

# End of edgesFromZScoreClustering
