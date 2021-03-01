'''
python3 script to generate a maxquant xml file from a template xml and a folder of raw files

Changes the following entries in the xml file:
    input file related: filePaths, experiments, fractions, ptms, paramGroupIndices, referenceChannel
    
Can also change proteome (--proteome): fastaFilePath and output folder (--output): fixedCombinedFolder

Author: Xuebing Wu (xw2629@columbia.edu)

TODO
    1. add option --threads
    2. add option to change sample name and fraction name

'''

import argparse
import os
from glob import glob
import sys

def valid_raw_folder(raw_dir):
    # check if the raw_dir folder exists and has *.raw files
    # return the number of *.raw files. -1: folder not exist
    nSample=-1
    if os.path.isdir(raw_dir):
        nSample=0
        files = os.listdir(raw_dir)
        for file in files:
            if file[-4:] == ".raw":
                nSample = nSample + 1
        if nSample == 0:
            print("WARNING: no *.raw files found in the input folder "+raw_dir,file=sys.stderr)
    else:
        print("WARNING: raw file folder doesn't exist: "+raw_dir,file=sys.stderr)
    return nSample
        
def generate_xml(template_xml,raw_dir,output_dir,proteome):
    nSample = valid_raw_folder(raw_dir)
    if nSample < 1: 
        return nSample # raw_dir not exist or has no *.raw files
    if not os.path.isfile(template_xml):
        return -2 # template xml file not found
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    fin=open(template_xml,'r')
    fout=open(os.path.abspath(output_dir)+'/mqpar.xml','w')
    skip=False
    for line in fin:
        if '<fastaFilePath>' in line:
            if proteome != 'NA':
                fout.write("<fastaFilePath>"+os.path.abspath(proteome)+"</fastaFilePath>\n")
            else:
                fout.write(line)
        elif '<fixedCombinedFolder>' in line:
            fout.write("<fixedCombinedFolder>" + os.path.abspath(output_dir) + "</fixedCombinedFolder>\n")
        elif '<filePaths>' in line:
            fout.write(line)
            skip=True
            for file in glob(raw_dir+'/*.raw'):
                fout.write("\t\t<string>"+os.path.abspath(file)+"</string>\n")
            fout.write("\t</filePaths>\n")
            fout.write("\t<experiments>\n")
            for i in range(nSample):
                fout.write("\t\t<string>Sample</string>\n")
            fout.write("\t</experiments>\n")
            fout.write("\t<fractions>\n")
            for i in range(nSample):
                #fout.write("\t\t<short>32767</short>\n")
                fout.write("\t\t<short>"+str(i+1)+"</short>\n")
            fout.write("\t</fractions>\n")
            fout.write("\t<ptms>\n")
            for i in range(nSample):
                fout.write("\t\t<boolean>False</boolean>\n")
            fout.write("\t</ptms>\n")
            fout.write("\t<paramGroupIndices>\n")
            for i in range(nSample):
                fout.write("\t\t<int>0</int>\n")
            fout.write("\t</paramGroupIndices>\n")
            fout.write("\t<referenceChannel>\n")
            for i in range(nSample):
                fout.write("\t\t<string></string>\n")
            fout.write("\t</referenceChannel>\n")
        elif '</referenceChannel>' in line:
            skip = False
        elif not skip:
            fout.write(line)
    fin.close()
    fout.close()
    print("- MaxQuant xml file created: "+os.path.abspath(output_dir)+'/mqpar.xml')
    return nSample

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--raw", dest='raw_dir', action='store',default='NA',
                         help='Rquired. Path to a folder with raw data files *.raw')

    parser.add_argument("--xml", dest='template_xml', action='store',default='NA',
                         help='Required. A template xml file')

    parser.add_argument("--output", dest='output_dir', action='store',default='NA',
                         help='Output directory. Default: same as input')

    parser.add_argument("--proteome", dest='proteome', action='store',default='NA',
                         help='Path to proteome fasta file. Default: use the one in the template xml')

    args = parser.parse_args()

    if args.raw_dir == 'NA' or args.template_xml =='NA':
        print("ERROR: You must specify path to the folder with the raw data (--raw) AND a template xml file (--xml)",file=sys.stderr)
        exit()
    
    if not os.path.isfile(args.template_xml):
        print("ERROR: the template xml file not found: "+args.template_xml,file=sys.stderr)
        exit()
        
    if args.output_dir == 'NA':
        args.output_dir = args.raw_dir
    
    generate_xml(args.template_xml,args.raw_dir,args.output_dir,args.proteome)