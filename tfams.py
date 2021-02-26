'''
A pipeline for detecting amino acid substitutions from proteomics data

Core algorithm/scripts from https://github.com/ernestmordret/substitutions/ by Ernest Mordret

Xuebing Wu

'''

# Required: please update the path to MaxQuant binary
MaxQuantCmd="dotnet /home/xw2629/software/MaxQuant/bin/MaxQuantCmd.exe"

# Option: setup default files
proteome='./reference/human.protein.fa'
transcriptome='./reference/human.CDS.fa'
template_xml='./template_xml/mqpar-human.xml'

import argparse
import os
from generate_xml import *

parser = argparse.ArgumentParser()

parser.add_argument("--input", dest='input_dir', action='store',default='NA',required=True,
                     help='Required. Path to a folder with raw data (*.raw)')

parser.add_argument("--output", dest='output_dir', action='store',default='NA',
                     help='Output folder name. Default: same as input')

parser.add_argument("--proteome", dest='proteome', action='store',default=proteome,
                     help='Path to proteome fasta file')

parser.add_argument("--transcriptome", dest='transcriptome', action='store',default=transcriptome,
                     help='Path to transcriptome (CDS only) fasta file')

parser.add_argument("--xml", dest='template_xml', action='store',default=template_xml,
                     help='A template xml file')
    
args = parser.parse_args()


if args.transcriptome == 'ecoli':
    args.transcriptome = os.path.abspath('./reference/ecoli.CDS.fa')

if args.proteome == 'ecoli':
    args.proteome = os.path.abspath('./reference/ecoli.protein.fa')

if args.input_dir == 'NA':
    print("ERROR: You must specify path to the folder with the raw data: --input")
    exit(1)
else:
    args.input_dir = os.path.abspath(args.input_dir)
    nSample = valid_raw_folder(args.input_dir) # will exit if no *.raw files found

if args.output_dir == 'NA':
    args.output_dir = args.input_dir
    
# create params.py
f = open("params.py", "w")
f.write("input_dir = '"+args.input_dir+"'\n")
f.write("output_dir = '"+args.output_dir+"'\n")
f.write("transcriptome = '"+args.transcriptome+"'\n")
f.close()    
os.system('cat params.core >> params.py')

print("Path to files:")
print("- proteome      : "+args.proteome)
print("- transcriptome : "+args.transcriptome)
print("- template xml  : "+args.template_xml)
print("- input         : "+args.input_dir)
print("- output        : "+args.output_dir)

if not os.path.isdir(args.output_dir):
    print('Creating output directory: '+args.output_dir)
    os.mkdir(args.output_dir)

if os.path.isfile(args.output_dir+'/combined/txt/allPeptides.txt'):
    print("MaxQuant search results detected in the output folder "+args.output_dir)
    print("- to run MaxQuant again, please change output directory")
else:
    print("Generate MaxQuant parameter file (xml)")
    generate_xml(args.template_xml,args.input_dir,args.output_dir,args.proteome)
    print("MaxQuant search of dependent peptides")
    cmd = MaxQuantCmd+ " " + args.output_dir + "/mqpar.xml"
    print('- '+cmd)
    os.system(cmd)
print("Detection and filtering")
import detect
import quantify
import plot