'''
A pipeline for detecting frameshift events from proteomics data

Xuebing Wu

'''

# Required: please update the path to MaxQuant binary
MaxQuantCmd="dotnet /home/xw2629/software/MaxQuant/bin/MaxQuantCmd.exe"

# Option: setup default files
transcriptome='./reference/human.CDS.fa'
template_xml='./template_xml/mqpar-frameshift.xml'

import argparse
import os
from generate_xml import *
from Bio import SeqIO
import time
import pandas as pd


def transcriptome2proteome(transcriptome):
    # note that this doesn't include 3' UTR
    if not os.path.isfile(transcriptome):
        return -1 # file not found
    out0 = open(transcriptome+'-proteome-frame0.fa','w')
    out1 = open(transcriptome+'-proteome-frame1.fa','w')
    out2 = open(transcriptome+'-proteome-frame2.fa','w')
    for record in SeqIO.parse(open(transcriptome,'rU'),'fasta'):
        record.seq = record.seq.upper()    
        if len(record.seq)%3 == 0 and record.seq[:3] in {'ATG','GTG','TTG','ATT','CTG'} and record.seq[-3:].translate()=='*':
            name=record.description.split(' ')[0]
            header = '>'+name+'|'+name+'|'+name
            translation = str(record.seq.translate()[:-1])#.replace('*','R')
            if len(translation)>0:
                out0.write(header+'\n'+translation+'\n')
            translation = str(record.seq[1:].translate()[:-1])#.replace('*','R')
            if len(translation)>0:
                out1.write(header+'\n'+translation+'\n')
            translation = str(record.seq[2:].translate()[:-1])#.replace('*','R')
            if len(translation)>0:
                out2.write(header+'\n'+translation+'\n')
    out0.close()
    out1.close()
    out2.close()
    return 0

def frameshift_detection(input_dir,output_dir,transcriptome,template_xml,frame0_peptides):
    if not ( os.path.isfile(transcriptome+'-proteome-frame0.fa') and os.path.isfile(transcriptome+'-proteome-frame1.fa') and os.path.isfile(transcriptome+'-proteome-frame2.fa')):
        print('Creating proteome from transcriptome: '+transcriptome)
        transcriptome2proteome(transcriptome)
    
    # skp frame0 run if path to peptides.txt is provided
    if os.path.isfile(frame0_peptides):
        n = 2
    else:
        n = 3
        
    for i in range(n):
        i=2-i # run frame 2, then frame 1, then frame 0 (slowest)
        if not os.path.isfile(output_dir+'/frame'+str(i)+'/combined/proc/Finish_writing_tables 11.finished.txt'):
            print("MaxQuant search in frame "+str(i))
            print("---------------------------------")
            if not os.path.isdir(output_dir+'/frame'+str(i)):
                os.mkdir(output_dir+'/frame'+str(i))
            generate_xml(template_xml,input_dir,output_dir+'/frame'+str(i),transcriptome+'-proteome-frame'+str(i)+'.fa')
            cmd = MaxQuantCmd+ " " + output_dir + "/frame"+str(i)+"/mqpar.xml "# "> "+args.output_dir + "/frame"+str(i)+"/log.txt 2>&1"
            print(cmd)
            os.system(cmd)
            os.system('wc -l '+output_dir+'/frame'+str(i)+'/combined/txt/peptides.txt')
        else:
            print("MaxQuant analysis has been completed for frame "+str(i))

    # load frame0 peptides
    if os.path.isfile(frame0_peptides):
        pep0 = pd.read_csv(frame0_peptides,sep = '\t',index_col='Sequence')
    else:
        pep0 = pd.read_csv(output_dir+'/frame0/combined/txt/peptides.txt',sep = '\t',index_col='Sequence')
        
    pep1 = pd.read_csv(output_dir+'/frame1/combined/txt/peptides.txt',sep = '\t',index_col='Sequence')
    pep2 = pd.read_csv(output_dir+'/frame2/combined/txt/peptides.txt',sep = '\t',index_col='Sequence')

    print('peptides identified in frame 0: '+str(len(pep0)))
    print('median intensity / 1e6: '+str(pep0['Intensity'].median()/1e6))

    print('peptides identified in frame 1: '+str(len(pep1)))

    allpep=' '.join(pep0.index)

    for seq in pep1.index:
        if seq in allpep:
            pep1=pep1.drop(seq)
    print('frameshift peptides in frame 1: '+str(len(pep1)))
    print('median intensity / 1e6: '+str(pep1['Intensity'].median()/1e6))

    print('peptides identified in frame 2: '+str(len(pep2)))

    for seq in pep2.index:
        if seq in allpep:
            pep2=pep2.drop(seq)

    print('frameshift peptides in frame 2: '+str(len(pep2)))
    print('median intensity / 1e6: '+str(pep2['Intensity'].median()/1e6))

    pep1.to_csv(output_dir+'/frame1/frameshift-peptides.txt')
    pep2.to_csv(output_dir+'/frame2/frameshift-peptides.txt')
                                                                                
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    parser.add_argument('input_dir',help='Required. Path to a folder with raw data (*.raw)')

    parser.add_argument("--output-dir",action='store',default='NA',
                         help='Output folder name. Default: same as input')

    parser.add_argument("--transcriptome", action='store',default=transcriptome,
                         help='Path to transcriptome (CDS only) fasta file')

    parser.add_argument("--template-xml", dest='template_xml', action='store',default=template_xml,
                         help='A template xml file')

    parser.add_argument("--frame0", dest='frame0', action='store',default="",
                         help='Path to existing peptides.txt from MaxQuant search in non-frameshifted proteome')

    args = parser.parse_args()


    if args.input_dir == 'NA':
        print("ERROR: You must specify path to the folder with the raw data: --input")
        exit(1)
    else:
        args.input_dir = os.path.abspath(args.input_dir)
        nSample = valid_raw_folder(args.input_dir)
        if nSample < 0:
            exit()

    if args.output_dir == 'NA':
        args.output_dir = args.input_dir

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    args.output_dir = args.output_dir +'/frameshift'
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    print("Path to files:")
    print("- transcriptome : "+os.path.abspath(args.transcriptome))
    print("- template xml  : "+os.path.abspath(args.template_xml))
    print("- output        : "+os.path.abspath(args.output_dir))
    print("- input         : "+os.path.abspath(args.input_dir))
    print("                  " + str(nSample) + " raw file(s)")
    print("- frame0        : "+os.path.abspath(args.frame0))

    frameshift_detection(args.input_dir,args.output_dir,args.transcriptome,args.template_xml,args.frame0)