#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:58:34 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import os
import argparse
from Bio.Align.Applications import MafftCommandline
from time import time

#%% Function definition
def align_mafft(BUSCODIR, OUTDIR, CONFIG_FILE=False):
	try:#check if input directory with fasta files exists
		files = os.listdir(BUSCODIR)
		print("%s busco fasta files found in buscodir." % str(len(files)))
	except:
		raise ValueError("buscodir %s can not be found." % BUSCODIR)
	try:#check if output directory exists and (if not) can be created
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
	except:
		raise ValueError("Outputdir can not be created in %s." % OUTDIR)

	print("Starting Mafft alignments...")
	for busco in files:
		infile = os.path.join(BUSCODIR, busco)
		settings = parse_maff_config(CONFIG_FILE)
		mafft_cmline = generate_cmdline(infile, settings)
		stdout, stderr = mafft_cmline()
		#outfile = os.path.join(OUTDIR, busco[:-4] + ".aln.fa")
		outfile = os.path.join(OUTDIR, busco.split('.')[0] + ".aln.fa")
		with open(outfile, "wt") as handle:
			handle.write(stdout)
#end
def generate_cmdline(infile, settings):
	mafft_cmline = MafftCommandline(settings['mafft_bin'], input=infile)
	#
	if settings['align_method'] == 'auto':
		#print("Running Mafft in auto mode.")
		mafft_cmline.set_parameter("--auto", True)
		mafft_cmline.set_parameter("--thread", 8)
		#
	else:
		if settings['align_method'] == 'AOM1':
			#print("Running Mafft in Accuracy-oriented mode: %s." % settings['align_method'])
			#AOM1: mafft --localpair --maxiterate 1000 input [> output]
			mafft_cmline.set_parameter("--localpair", True)
			mafft_cmline.set_parameter("--maxiterate", 1000) #mafft_call.maxiterate = 1000
			#
		elif settings['align_method'] == 'AOM2':
			#print("Running Mafft in Accuracy-oriented mode: %s." % settings['align_method'])
			#AOM2: mafft --globalpair --maxiterate 1000 input [> output]
			mafft_cmline.set_parameter("--globalpair", True)
			mafft_cmline.set_parameter("--maxiterate", 1000) #mafft_call.maxiterate = 1000
			#
		elif settings['align_method'] == 'AOM3':
			#print("Running Mafft in Accuracy-oriented mode: %s." % settings['align_method'])
			#AOM3: mafft --ep 0 --genafpair --maxiterate 1000 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 1000)
			mafft_cmline.set_parameter("--genafpair", True)
			mafft_cmline.set_parameter("--ep", 0)
			#
		elif settings['align_method'] == 'SOM1':
			#print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM1: mafft --retree 2 --maxiterate 2 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 2)
			mafft_cmline.set_parameter("--retree", 2)
			#
		elif settings['align_method'] == 'SOM2':
			#print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM2: mafft --retree 2 --maxiterate 1000 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 2)
			mafft_cmline.set_parameter("--retree", 2)
			#
		elif settings['align_method'] == 'SOM3':
			#print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM3: mafft --retree 2 --maxiterate 0 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 0)
			mafft_cmline.set_parameter("--retree", 2)
			#
		elif settings['align_method'] == 'SOM4':
			#print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM4: mafft --retree 1 --maxiterate 0 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 0)
			mafft_cmline.set_parameter("--retree", 1)
			#
		elif settings['align_method'] == 'SOM5':
			#print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM5: mafft --retree 2 --maxiterate 2 --nofft input [> output]
			mafft_cmline.set_parameter("--maxiterate", 2)
			mafft_cmline.set_parameter("--retree", 2)
			mafft_cmline.set_parameter("--nofft", True)
			#
		elif settings['align_method'] == 'SOM6':
			#print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM6: mafft --retree 2 --maxiterate 0 --nofft input [> output]
			mafft_cmline.set_parameter("--maxiterate", 0)
			mafft_cmline.set_parameter("--retree", 2)
			mafft_cmline.set_parameter("--nofft", True)
			#
		elif settings['align_method'] == 'SOM7':
			#print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM7: mafft --retree 1 --maxiterate 0 --nofft --parttree input [> output]
			mafft_cmline.set_parameter("--maxiterate", 0)
			mafft_cmline.set_parameter("--retree", 1)
			mafft_cmline.set_parameter("--nofft", True)
			mafft_cmline.set_parameter("--parttree", True)
			#
		elif settings['align_method'] == 'manual':
			#print("Running Mafft in manual mode.")
			mafft_cmline.set_parameter("--maxiterate", int(settings['maxiterate']))
			if settings['type_method'] == 'accuracy':
				if settings['pairwise_method'] == 'genafpair':
					mafft_cmline.set_parameter("--genafpair", True)
					mafft_cmline.set_parameter("--ep", int(settings['ep']))
				elif settings['pairwise_method'] == 'globalpair':
					mafft_cmline.set_parameter("--globalpair", True)
				elif settings['pairwise_method'] == 'localpair':
					mafft_cmline.set_parameter("--localpair", True)
			elif settings['type_method'] == 'speed':
				mafft_cmline.set_parameter("--retree", int(settings['retree']))
				mafft_cmline.set_parameter("--nofft", int(settings['nofft']))
				mafft_cmline.set_parameter("--parttree", int(settings['parttree']))
			#
		if settings['output_setting'] == 'manual':
			#print("entró en output setting manual")
			# Output options
			mafft_cmline.set_parameter("--thread", int(settings['threads']))
			mafft_cmline.set_parameter("--clustalout", int(settings['clustalout']))
			mafft_cmline.set_parameter("--inputorder", int(settings['inputorder']))
			mafft_cmline.set_parameter("--reorder", int(settings['reorder']))
			mafft_cmline.set_parameter("--treeout", int(settings['treeout']))
			mafft_cmline.set_parameter("--quiet", int(settings['quiet']))
	return mafft_cmline
#end
def parse_maff_config(CONFIG_FILE):
	settings = {}
	if not CONFIG_FILE:
		settings['align_method'] = 'auto'
		settings['mafft_bin'] = 'mafft'
	else:
		try:
			with open(CONFIG_FILE, "rt") as INPUT:
				for line in INPUT:
					if line[0] == '#' or line[0] == '\n' or line[0] == ' ':
						continue
					line = line.strip()
					param,value = line.split('=')
					settings[param.rstrip()] = value.strip().split(" ")[0]
		except:
				raise ValueError("Config file %s can not be found." % CONFIG_FILE)
	return settings
#end
def check_arguments(BUSCODIR, OUTDIR, CONFIG_FILE=False):
	if not os.path.isdir(BUSCODIR):
		raise ValueError("buscodir with fasta files %s does not exists." % BUSCODIR)
	try:#check if output directory exists and (if not) can be created
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
	except:
		raise ValueError("Outputdir can not be created in %s." % OUTDIR)
	if CONFIG_FILE:
		if not os.path.isfile(CONFIG_FILE):
			raise ValueError("Config file %s does not exists." % CONFIG_FILE)
#end
#%% Menu -> is executed when the script is called independently
def usage():
	parser = argparse.ArgumentParser(
		description='''02_align_concatenate.py takes FASTA files with BUSCO orthologs, generate the alignment of each one and concatenate all alignments into a matrix.''',
		epilog="""End of the help""")
	parser.add_argument('-b', '--buscodir', type=str, required=True, help='Path to the directory containing BUSCO FASTA files with the orthologous sequences of all genomes, e.g. "common_busco_seqs" directory part of 01_find_singlecopy.py output.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Full or relative path to the directory where the aligments, phylogenetic matrix and coordinates will be placed. If it does not exists, it will be created.')
	parser.add_argument('-c', '--config', metavar='<config file>', type=str, required=False, help='Config file for alignment setting. You can find a template config file in the docs directory. If no file is provided, the alignments will be done using default parameters.')
	#parser.add_argument('-a', '--aligner', type=str, required=False , choices=['mafft', 'muscle', 'clustalw'], default='mafft', help='Alignment program to use.')
	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
	args = usage()
	print("Checking paths...")
	start = time()
	check_arguments(args.buscodir, args.outdir, args.config)
	print("Starting Mafft alignments...")
	align_mafft(args.buscodir, args.outdir, args.config)
	print("Alignments completed...")
	print(f'Time taken to run: {time() - start} seconds.')
#end







