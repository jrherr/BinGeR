BinGeR 
version: 0.0.1
author: Chengwei Luo

BinGeR: in-situ Genome recovery tool for series metagenomes
==========================

BinGeR: in-situ Genome recovery tool for series metagenomes

BinGeR (Binner for Genome Recovery) is a platform for de novo genome recovery from series metagenomes. It integrates information from single copy gene content, contig coverage correlation, and contig oligo-nucleotide composition correlation, and contig alignments to bin genome fragments together. It is written in Python, therefore it should run on Windows, Mac OS, and Unix/Linux.

Reference: [manuscript in preparation]

Copyright: Chengwei Luo, Konstantinidis Lab, Georgia Institute of Technology, 2013

Contact: Chengwei Luo (Email: luo.chengwei@gatech.edu)

[How to install]
==========================
Dependencies:

Python >= 2.6
BLAT 
prodigal
hmmscan

Python libraries:

NetworkX >= 1.7 
Pysam >= 0.6
Numpy >= 1.5.1
Scipy >= 0.12.0
BioPython >= 1.58
scikit-learn >= 0.14

Packages with older versions might work, but not tested. You can certainly manually install all of the packages, however, I recommend using Anaconda, which is a Python distribution for big data processing, it bundles many packages for data scientists. For more information, go to:

https://store.continuum.io/cshop/anaconda/

Anaconda includes all the dependencies to run BinGeR.

To install, unzip/untar/uncompress the BinGeR package, and in terminal cd to the unpacked directory, and run:

$ python setup.py install

You can alway supply --prefix to alter the installation directory:

$ python setup.py install --prefix = user/defined/installation/directory

If you have everything ready, this should take care of the installation.

[How to start a BinGeR project]
==================================
Here is what you need to get started:

1, a text file listing the sample names, one per line;

2, pre-contigs for each sample in fasta format. By pre-contigs, I mean assemblies you manage to get from standard assemblers such as Velvet, SOAPdenovo, etc. By default, these files are stored at 'Assemblies/', otherwise when running BinGeR, you'll need to use -a/--assemblies_dir to specify where they are;

3, all-vs-all reads mapping to these pre-contigs in BAM format, sorted and indexed. For instance, if you have 4 samples, you will need 4x4 BAM files. If you are not familiar with how to generate BAM files, refer to Samtools page: http://samtools.sourceforge.net.

With all these, you need to do the following using the scripts in the 'utils' directory (they are easy to be run in a distributed manner, so if you have access to a cluster, this will save you time):

1, generate coverage files using utils/bamCoverage.py; by default they are saved in the 'Coverage/' directory;

2, generate 3-mer and 4-mer stats using utils/oligoNtZScore.py; by default they are saved in the 'ZScores/' directory;

3, generate hmmscan files using utils/runHMMScan.py; by default they are saved in the 'HMMScan/' directory;

4, save all the 'on-diagonal' BAM files into a folder (.bam and .bai), which by default is 'Bams/', however, you can always specify the location of the folder when running BinGeR otherwise. (you don't need the "off-diagonal" BAM files from now on). By "on-diagonal" I mean self-vs-self mappings, for example 'sampleA.vs.sampleA.sorted.bam'.

Now you are ready to run BinGeR.py

[How to run BinGeR]
==========================
After running setup.py, you may want to either open a new tab in the terminal, or in current terminal tab:
$ source ~/.profile

In terminal, as simple as:
$ BinGeR.py -l <sample.txt> -o <output_dir> [options]

For more detailed settings, please type:

$ BinGeR.py --help

This should print out a detailed settings for you.

Note: to start with a project, you will need to use the scripts in the utils directory to produce the files first


[Contents of db folder]
==========================
In db/, you will find a few db files:
HMM.txt
ncbiNodes.lib
ncbiSciNames.lib
singleCopy.prot.tar.gz

[Contents of utils folder]
==========================
In utils/, you will see a few Python scripts, as listed below:

bamCoverage.py
oligoNtZScore.py
runHMMScan.py

They are utility scripts that will assist you in preparing the input files necessary for BinGeR.

For detailed usage of each utils script, please run:

$ python utils/script.py --help

There is also an hmm file, "SingleCopyGenes.HMM", which is used as the database for runHMMScan.py
