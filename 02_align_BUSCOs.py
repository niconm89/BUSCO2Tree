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
import subprocess

#%% Function definition
def align_config_mafft(BUSCODIR, OUTDIR, CONFIG_FILE=False):
	'Receive directory with common BUSCO multi FASTA files, output directory and config file (optional) and generate alignments using MAFFT.'
	try: #check if input directory with fasta files exists
		files = os.listdir(BUSCODIR)
		print("%s busco fasta files found in buscodir." % str(len(files)))
	except:
		raise ValueError("buscodir %s can not be found." % BUSCODIR)
	try: #check if output directory exists, otherwise it will be created
		if not os.path.isdir(OUTDIR):
			os.mkdir(OUTDIR)
	except:
		raise ValueError("Outputdir can not be created in %s." % OUTDIR)

	for busco in files:
		infile = os.path.join(BUSCODIR, busco) #multi FASTA file to align
		settings = parse_mafft_config(CONFIG_FILE) #reading alignment settings
		mafft_cmline = generate_cmdline(infile, settings) #generating MAFFT command
		stdout, stderr = mafft_cmline() #running MAFFT
		outfile = os.path.join(OUTDIR, busco.split('.')[0] + ".aln.fa") #naming alignment file
		with open(outfile, "wt") as handle: #saving alignment file
			handle.write(stdout)
#end
def generate_cmdline(infile, settings):
	mafft_cmline = MafftCommandline(settings['mafft_bin'], input=infile)
	#
	if settings['align_method'] == 'auto':
		print("Running Mafft in auto mode.")
		mafft_cmline.set_parameter("--auto", True)
		mafft_cmline.set_parameter("--thread", 8)
		#
	else:
		if settings['align_method'] == 'AOM1':
			print("Running Mafft in Accuracy-oriented mode: %s." % settings['align_method'])
			#AOM1: mafft --localpair --maxiterate 1000 input [> output]
			mafft_cmline.set_parameter("--localpair", True)
			mafft_cmline.set_parameter("--maxiterate", 1000) #mafft_call.maxiterate = 1000
			#
		elif settings['align_method'] == 'AOM2':
			print("Running Mafft in Accuracy-oriented mode: %s." % settings['align_method'])
			#AOM2: mafft --globalpair --maxiterate 1000 input [> output]
			mafft_cmline.set_parameter("--globalpair", True)
			mafft_cmline.set_parameter("--maxiterate", 1000) #mafft_call.maxiterate = 1000
			#
		elif settings['align_method'] == 'AOM3':
			print("Running Mafft in Accuracy-oriented mode: %s." % settings['align_method'])
			#AOM3: mafft --ep 0 --genafpair --maxiterate 1000 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 1000)
			mafft_cmline.set_parameter("--genafpair", True)
			mafft_cmline.set_parameter("--ep", 0)
			#
		elif settings['align_method'] == 'SOM1':
			print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM1: mafft --retree 2 --maxiterate 2 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 2)
			mafft_cmline.set_parameter("--retree", 2)
			#
		elif settings['align_method'] == 'SOM2':
			print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM2: mafft --retree 2 --maxiterate 1000 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 2)
			mafft_cmline.set_parameter("--retree", 2)
			#
		elif settings['align_method'] == 'SOM3':
			print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM3: mafft --retree 2 --maxiterate 0 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 0)
			mafft_cmline.set_parameter("--retree", 2)
			#
		elif settings['align_method'] == 'SOM4':
			print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM4: mafft --retree 1 --maxiterate 0 input [> output]
			mafft_cmline.set_parameter("--maxiterate", 0)
			mafft_cmline.set_parameter("--retree", 1)
			#
		elif settings['align_method'] == 'SOM5':
			print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM5: mafft --retree 2 --maxiterate 2 --nofft input [> output]
			mafft_cmline.set_parameter("--maxiterate", 2)
			mafft_cmline.set_parameter("--retree", 2)
			mafft_cmline.set_parameter("--nofft", True)
			#
		elif settings['align_method'] == 'SOM6':
			print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM6: mafft --retree 2 --maxiterate 0 --nofft input [> output]
			mafft_cmline.set_parameter("--maxiterate", 0)
			mafft_cmline.set_parameter("--retree", 2)
			mafft_cmline.set_parameter("--nofft", True)
			#
		elif settings['align_method'] == 'SOM7':
			print("Running Mafft in Speed-oriented mode: %s " % settings['align_method'])
			#SOM7: mafft --retree 1 --maxiterate 0 --nofft --parttree input [> output]
			mafft_cmline.set_parameter("--maxiterate", 0)
			mafft_cmline.set_parameter("--retree", 1)
			mafft_cmline.set_parameter("--nofft", True)
			mafft_cmline.set_parameter("--parttree", True)
			#
		elif settings['align_method'] == 'manual':
			print("Running Mafft in manual mode.")
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
			# Output options in manual mode
			mafft_cmline.set_parameter("--thread", int(settings['threads']))
			mafft_cmline.set_parameter("--clustalout", int(settings['clustalout']))
			mafft_cmline.set_parameter("--inputorder", int(settings['inputorder']))
			mafft_cmline.set_parameter("--reorder", int(settings['reorder']))
			mafft_cmline.set_parameter("--treeout", int(settings['treeout']))
			mafft_cmline.set_parameter("--quiet", int(settings['quiet']))
	return mafft_cmline
#end
def parse_mafft_config(CONFIG_FILE):
	'Parse the MAFFT config file.'
	settings = {} #dictionary to save all settings
	if not CONFIG_FILE: #if no config file is provided, MAFFT will be run in auto mode.
		settings['align_method'] = 'auto'
		settings['mafft_bin'] = 'mafft' #MAFFT is assumed to be available in the $PATH variable.
	else: #config file is passed
		try:
			with open(CONFIG_FILE, "rt") as INPUT:
				for line in INPUT:
					if line[0] == '#' or line[0] == '\n' or line[0] == ' ': #avoiding lines
						continue
					line = line.strip()
					param,value = line.split('=')
					settings[param.rstrip()] = value.strip().split(" ")[0]
		except:
				raise ValueError("Config file %s can not be found." % CONFIG_FILE)
	return settings
#end
def align_command_mafft(BUSCODIR, OUTDIR, COMMAND):
    'Receive directory with common BUSCO multi FASTA files, output directory and a MAFFT command (file) to generate alignments.' 
    if not COMMAND:
        raise ValueError("Command for MAFFT alignment is empty.")
    base_cmd_mafft = COMMAND
    for busco in os.listdir(BUSCODIR):
        infile = os.path.join(BUSCODIR, busco)
        outfile = os.path.join(OUTDIR, busco.split(".")[0] + ".mafft.fa")
        cmd_mafft = base_cmd_mafft + " " + infile + " > " + outfile
        #print(cmd_mafft)
        #subprocess.call([cmd_mafft], shell=True)
        run_mafft = subprocess.call([cmd_mafft], shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
#end
def trim_alns(OUTDIR, TRIMPARAMS):
    'Receive parameters to run TrimAl, avoiding input and output files.'
    #do something
    list_of_alns = os.listdir(OUTDIR)
    for aln in list_of_alns:
        infile = os.path.join(OUTDIR, aln)
        outfile = os.path.join(OUTDIR, aln.split(".")[0] + ".aln.trimmed.fa")
        cmd_trimal = "trimal -in " + infile + " -out " + outfile + " "
        if TRIMPARAMS:
            cmd_trimal += TRIMPARAMS
            subprocess.call([cmd_trimal], shell=True)
#end
def check_arguments(BUSCODIR, OUTDIR, CONFIG_FILE=False, COMMAND_FILE=False):
	'Checks argumments of 02_align_BUSCOs.py script.'
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
		description='''02_align_concatenate.py takes multi FASTA files with common BUSCO groups (orthologs), generate the alignment of each one usign MAFFT and concatenate all alignments into a matrix.''',
		epilog="""End of the help""")
	parser.add_argument('-b', '--buscodir', type=str, required=True, help='Path to the directory containing BUSCO multi FASTA files with the ortholog sequences of all genomes, e.g. "common_busco_seqs" directory that is created as part of the 01_find_singlecopy.py output.')
	parser.add_argument('-o', '--outdir', type=str, required=True, help='Path to the directory where the aligments, phylogenetic matrix and coordinates will be placed. If it does not exists, it will be created.')
	parser.add_argument('-c', '--config', metavar='<config file>', type=str, required=False, help='Configuration file for alignment settings. You can find a template of this file in the example directory. If no file is provided, the alignments will be done using default parameters of MAFFT.')
	parser.add_argument('-m', '--command', metavar='<command file>', type=str, required=False, help='Command (between quote marks) to use in the MAFFT alingment. The names of input and out files must be avoided, e.g. "mafft --thread 8 --unalignlevel 0.1 --leavegappyregion --ep 0.12 --globalpair --maxiterate 1000".')
	parser.add_argument('-t', '--trim', action='store_true', required=False, help='Trim alignments using TrimAl in automated mode, which must be available in the path.')
	parser.add_argument('-p', '--trimparams', required=False, help='TrimAl parameters to use. They must be passed in quotes marks and avoiding input and output files, e.g. "-gt 0.3 -nogaps -phylip".')

	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
	args = usage() #reading arguments
	start = time() #time 0
	print("Step 0: Checking paths...")
	check_arguments(args.buscodir, args.outdir, args.config, args.command)
	print("Step 1: Starting Mafft alignments...")
	if args.command:
		print("MAFFT will be excecuted using a command introduced by the user.")
		align_command_mafft(args.buscodir, args.outdir, args.command)
	else:
		align_config_mafft(args.buscodir, args.outdir, args.config)
	print("Alignments completed...")
	if args.trim: #trimming alingments
		print("Step 2: Trimming alignments...")
		trim_alns(args.outdir, args.trimparams)
	print(f'Time taken to run: {time() - start} seconds.')
#end