#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 13:33:47 2022

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import os

#%% Function definition
def dirpath(DIRPATH):
	if not os.path.isdir(DIRPATH):
		print("Error: BUSCO output directory does not exist. Check your paths! " + DIRPATH)
		exit(1)

def filepath(FILEPATH):
	if not os.path.isdir(FILEPATH):
		print("Error: BUSCO output directory does not exist. Check your paths! " + FILEPATH)
		exit(1)

def parentdir(PATH):
	if '/' in PATH:
		elements_path = PATH.split('/')
		parentdir = '/'.join(elements_path[:-1])
		if not os.path.isdir(parentdir):
			print("Error: Output directory can not be found or created in the path you have passed. Parent directory does not exist! " + parentdir)
			exit(1)
#end