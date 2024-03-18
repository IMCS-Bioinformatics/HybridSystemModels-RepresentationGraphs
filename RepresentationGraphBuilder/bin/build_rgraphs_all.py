################################################################################
# build_rgraphs_all.py
# Python script for building representation graphs for all HSM model state spaces in given directory
# Usage: 
# build_rgraphs_all.py <source dir> <destination dir> [-draw] [-ini]
# Requirements:
# Python 3.8.1 or newer
# Dependencies:
# numpy-1.24.4, netwrokx-3.1, igraph-0.11.4, N2G-0.3.3
# Contributors: 
# Institute of Mathematics and Computer Science, University of Latvia
# v_1.0.17, 18.03.2024
# Distributed under GPLv3 license
# Copyright (c) 2024 Juris Viksna
################################################################################

import os
import sys
import textwrap
import argparse
from argparse import ArgumentParser, HelpFormatter

################################################################################
# check and assign the input arguments
################################################################################

class RawFormatter(HelpFormatter):
    def _fill_text(self, text, width, indent):
        return "\n".join([textwrap.fill(line, width) for line in textwrap.indent(textwrap.dedent(text), indent).splitlines()])

parser = argparse.ArgumentParser(description="Computes representation graphs for all HSM model state spaces in given directory\n\nNote:\nAll files with '.txt' extension will be assumed to be HSM state space files\nEnsure that both input and output directories exist",formatter_class=RawFormatter)
parser.add_argument("-i",dest="in",help="input directory",metavar="<input dir>",required=True)
parser.add_argument("-o",dest="out",help="output directory",metavar="<output dir>",required=True)
parser.add_argument("-draw",dest="draw",action="store_true",help="draw also graphml files",required=False)
parser.add_argument("-ini",dest="ini",action="store_true",help="analyse only parts reachable from INI states",required=False)
args = vars(parser.parse_args())

indir = args['in']
outdir = args['out']

if not os.path.isdir(indir):
	print("Input directory not found!")
	exit()
if not os.path.isdir(outdir):
	print("Output directory not found!")
	exit()

	
all_files = os.listdir(indir)
hsm_files = [file for file in all_files if os.path.splitext(file)[1] == '.txt']
for infile in hsm_files:
	f_pref = infile.rsplit('.', maxsplit=1)[0]
	f_outgraph = f_pref+"_rg.txt"
	f_outdraw = f_pref+"_rg.graphml"
	pass_args = " -s -i "+indir+"/"+infile+" -o "+outdir+"/"+f_outgraph
	if args['draw']:
		pass_args = pass_args+" -draw "+outdir+"/"+f_outdraw
	if args['ini']:
		pass_args = pass_args+" -ini "
	print(pass_args)
	os.system("python build_rgraphs_file.py "+pass_args)

################################################################################
################################################################################
################################################################################