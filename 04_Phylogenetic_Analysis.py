#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 00:37:16 2022

@author: root
"""
#%% Imports
import subprocess

#%% Functionts definition

def model_partitions(MATRIXFILE, PARTITIONFILE, PREFIX, BOOTSTRAP, THREADS):
	#iqtree -s supermatrix.aln.faa.phy -p partitions-scheme.txt -m MFP 
	#--seqtype AA --prefix Drosophila -B 1000 --mem 50 -T 8
	
	IQTree_MFP = "/home/nmoreyra/Soft/miniconda3/envs/spyder/bin/iqtree -s " + MATRIXFILE + " -p " + PARTITIONFILE + " -m MFP --seqtype AA --prefix " + PREFIX + " -B " + str(BOOTSTRAP) + " -T " + str(THREADS)
	process = subprocess.Popen(IQTree_MFP.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()
#end