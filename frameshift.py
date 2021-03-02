'''
A pipeline for detecting frameshift events from proteomics data

Xuebing Wu

'''

# Required: please update the path to MaxQuant binary
MaxQuantCmd="dotnet /home/xw2629/software/MaxQuant/bin/MaxQuantCmd.exe"

# Option: setup default files
transcriptome='./reference/human.CDS.fa'
proteome='./reference/human.protein.fa'
template_xml='./template_xml/mqpar-frameshift.xml'

import argparse
import os
from generate_xml import *
from Bio import SeqIO
import time
import pandas as pd

def UTR_proteome(transcriptome):
    '''
    input: gencode coding gene sequence with UTRs
    '''
    if not os.path.isfile(transcriptome):
        return -1 # file not found
    out = open(transcriptome+'-UTR-proteome.fa','w')
    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
        record.seq = record.seq.upper()    
        header=record.description.split('|') #  >ENST00000641515.2|ENSG00000186092.6|OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-202|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|
        txpt_id = header[0]
        gene_symbol = header[5]
        txpt_len = int(header[6])
        if header[7][:3] == 'CDS':
            s,e = header[7][4:].split('-')
        elif header[8][:3] == 'CDS':
            s,e = header[8][4:].split('-')
        elif header[9][:3] == 'CDS':
            s,e = header[9][4:].split('-')
        else:
            pass
        s=int(s)
        e=int(e)
        # 5' UTR
        if s > 6:
            utr5seq = record.seq[:(s+18)]
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-0\n'+str(utr5seq[0:])+'\n')
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-1\n'+str(utr5seq[1:])+'\n')
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-2\n'+str(utr5seq[2:])+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-0\n'+str(utr5seq[0:].translate()[:-1])+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-1\n'+str(utr5seq[1:].translate()[:-1])+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-2\n'+str(utr5seq[2:].translate()[:-1])+'\n')
        if txpt_len - e > 6: # 3' UTR
            utr3seq = record.seq[(e-18):]
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-0\n'+str(utr3seq[0:])+'\n')
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-1\n'+str(utr3seq[1:])+'\n')
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-2\n'+str(utr3seq[2:])+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR3-0\n'+str(utr3seq[0:].translate()[:-1])+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR3-1\n'+str(utr3seq[1:].translate()[:-1])+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR3-2\n'+str(utr3seq[2:].translate()[:-1])+'\n')
    out.close()
    return 0
    

def transcriptome2proteome(transcriptome): # TODO: combine the two frame in a single reference
    '''
    generate protein sequences in three reading frames
    - input  : a fasta file containing CDS of all transcripts, each starts with start codon and ends with stop codon
    - output : three fasta files for proteins encoded in each reading frame. 
               - output in the same folder as input
               - output file suffix: -proteome-frame[0/1/2].fa
    '''
    # note that this doesn't include 3' UTR
    if not os.path.isfile(transcriptome):
        return -1 # file not found
    out0 = open(transcriptome+'-proteome-frame0.fa','w')
    out1 = open(transcriptome+'-proteome-frame1.fa','w')
    out2 = open(transcriptome+'-proteome-frame2.fa','w')
    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
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

def frameshift_detection(input_dir,output_dir,transcriptome,template_xml):
    if not ( os.path.isfile(transcriptome+'-proteome-frame0.fa') and os.path.isfile(transcriptome+'-proteome-frame1.fa') and os.path.isfile(transcriptome+'-proteome-frame2.fa')):
        print('Creating proteome from transcriptome: '+transcriptome)
        transcriptome2proteome(transcriptome)
        
    for i in [1,2]:
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
    

    pep1 = pd.read_csv(output_dir+'/frame1/combined/txt/peptides.txt',sep = '\t',index_col='Sequence')
    print('peptides identified in frame 1: '+str(len(pep1)))

    pep2 = pd.read_csv(output_dir+'/frame2/combined/txt/peptides.txt',sep = '\t',index_col='Sequence')
    print('peptides identified in frame 2: '+str(len(pep2)))

    print('load the proteome (frame 0)')
    allpep = open(transcriptome+'-proteome-frame0.fa').read()

    for seq in pep1.index:
        if seq in allpep:
            pep1=pep1.drop(seq)
    print('true frameshift peptides in frame 1 but not frame 0: '+str(len(pep1)))
    print('median intensity / 1e6: '+str(pep1['Intensity'].median()/1e6))


    for seq in pep2.index:
        if seq in allpep:
            pep2=pep2.drop(seq)

    print('true frameshift peptides in frame 2 but not frame 0: '+str(len(pep2)))
    print('median intensity / 1e6: '+str(pep2['Intensity'].median()/1e6))

    pep1.to_csv(output_dir+'/frame1/frameshift-peptides.txt')
    pep2.to_csv(output_dir+'/frame2/frameshift-peptides.txt')
                                                                                
def UTR_translation_detection(input_dir,output_dir,transcriptome,proteome,template_xml):
    if not os.path.isfile(transcriptome+'-UTR-proteome.fa'):
        print('Creating UTR proteome from transcriptome: '+transcriptome)
        UTR_proteome(transcriptome)
    else:
        print('Using existing UTR proteome: '+transcriptome+'-UTR-proteome.fa')
        
    if not os.path.isfile(output_dir+'/combined/txt/peptides.txt'):
        print("MaxQuant search for UTR peptides")
        print("---------------------------------")
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        generate_xml(template_xml,input_dir,output_dir,transcriptome+'-UTR-proteome.fa')
        cmd = MaxQuantCmd+ " " + output_dir + "/mqpar.xml "
        print(cmd)
        os.system(cmd)
        os.system('wc -l '+output_dir+'/combined/txt/peptides.txt')
    else:
        print("Found existing MaxQuant search result in UTRs: "+output_dir+'/combined/txt/peptides.txt')
    
    print("Removing candidate UTR peptides also found in the normal proteome")
    npep = filter_with_proteome(output_dir+'/combined/txt/peptides.txt',proteome)
    print(str(npep)+' UTR peptides saved to '+output_dir+'/combined/txt/peptides.txt-filtered.txt')
                  
def filter_with_proteome(path_to_pep,path_to_proteome):
    '''
    path_to_pep: path to peptides.txt
    path_to_proteome: path to reference proteome file
    '''
    pep = pd.read_csv(path_to_pep,sep = '\t',index_col='Sequence')
    allpep = open(path_to_proteome).read()
    for seq in pep.index:
        if seq in allpep:
            pep=pep.drop(seq)
    pep.to_csv(path_to_pep+'-filtered.txt')
    return len(pep)
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    parser.add_argument('input_dir',help='Required. Path to a folder with raw data (*.raw)')

    parser.add_argument("--output-dir",action='store',default='NA',
                         help='Output folder name. Default: same as input')

    parser.add_argument("--transcriptome", action='store',default=transcriptome,
                         help='Path to transcriptome (CDS only) fasta file')
              
    parser.add_argument("--proteome", action='store',default=proteome,
                         help='Path to proteome fasta file')

    parser.add_argument("--template-xml", dest='template_xml', action='store',default=template_xml,
                         help='A template xml file')

    parser.add_argument("--utr", action="store_true", 
                    help="run UTR peptide analysis instead of frameshift")
              
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

    if args.utr:
        args.output_dir = args.output_dir +'/utr'
    else:
        args.output_dir = args.output_dir +'/frameshift'
              
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    print("Path to files:")
    print("- transcriptome : "+os.path.abspath(args.transcriptome))
    print("- proteome : "+os.path.abspath(args.proteome))
    print("- template xml  : "+os.path.abspath(args.template_xml))
    print("- output        : "+os.path.abspath(args.output_dir))
    print("- input         : "+os.path.abspath(args.input_dir))
    print("                  " + str(nSample) + " raw file(s)")

    if args.utr:
        UTR_translation_detection(args.input_dir,args.output_dir,args.transcriptome,args.proteome,args.template_xml)
    else:
        frameshift_detection(args.input_dir,args.output_dir,args.transcriptome,args.template_xml)
