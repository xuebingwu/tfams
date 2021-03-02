'''
A pipeline for detecting frameshift events from proteomics data

Xuebing Wu

'''

# Required: please update the path to MaxQuant binary
MaxQuantCmd="dotnet /home/xw2629/software/MaxQuant/bin/MaxQuantCmd.exe"

# Option: setup default files
transcriptome_frameshift='./reference/human.CDS.fa'    # for frameshift
transcriptome_lncrna='./reference/human.lncRNA.fa'     # for lncRNA
transcriptome_mrna='./reference/human.mRNA.fa'         # for UTR
transcriptome_intron='./reference/human.intron.fa'     # for intron. downloaded from UCSC, gencode.v32, +9nt flanking sequence
proteome='./reference/human.protein.fa'
template_xml='./template_xml/mqpar-standard.xml'

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
    out = open(transcriptome+'-utr-proteome.fa','w')
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
        if s > 9: # at least 9nt 5' UTR
            utr5seq = record.seq[:(s+9)]
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-0\n'+str(utr5seq[0:])+'\n')
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-1\n'+str(utr5seq[1:])+'\n')
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-2\n'+str(utr5seq[2:])+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-0\n'+str(utr5seq[0:].translate())+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-1\n'+str(utr5seq[1:].translate())+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-2\n'+str(utr5seq[2:].translate())+'\n')
        if txpt_len - e >= 9: # 3' UTR at least 9 nt
            utr3seq = record.seq[(e-18):]
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-0\n'+str(utr3seq[0:])+'\n')
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-1\n'+str(utr3seq[1:])+'\n')
            #out.write('>'+txpt_id+'|'+gene_symbol+'|UTR5-2\n'+str(utr3seq[2:])+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR3-0\n'+str(utr3seq[0:].translate())+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR3-1\n'+str(utr3seq[1:].translate())+'\n')
            out.write('>'+txpt_id+'|'+gene_symbol+'|UTR3-2\n'+str(utr3seq[2:].translate())+'\n')
    out.close()
    return 0
    
def lncRNA_or_intron_proteome(transcriptome,noncoding_type): 
    '''
    generate protein sequences in three reading frames in lncRNAs
    - input  : a fasta file containing sequences of all lncRNAs
    - output : fasta files for proteins encoded in each reading frame, regardless of AUG
               - output in the same folder as input
               - output file suffix: -proteome.fa
    '''
    if not os.path.isfile(transcriptome):
        return -1 # file not found
    noncoding_type = str(noncoding_type).lower()
    
    out = open(transcriptome+'-'+noncoding_type+'-proteome.fa','w')

    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
        record.seq = record.seq.upper()    
        name=record.description
        out.write('>'+record.description+noncoding_type+'-frame-0\n'+str(record.seq[0:].translate())+'\n')
        out.write('>'+record.description+noncoding_type+'-frame-1\n'+str(record.seq[1:].translate())+'\n')
        out.write('>'+record.description+noncoding_type+'-frame-2\n'+str(record.seq[2:].translate())+'\n')
    out.close()
    return 0

def frameshift_proteome(transcriptome):
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
    out = open(transcriptome+'-frameshift-proteome.fa','w')
    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
        record.seq = record.seq.upper()    
        if len(record.seq)%3 == 0 and record.seq[:3] in {'ATG','GTG','TTG','ATT','CTG'} and record.seq[-3:].translate()=='*':
            description = record.description.split(' ')
            out.write('>'+description[0]+'|'+description[6]+'|frameshift+1\n'+str(record.seq[1:].translate())+'\n')
            out.write('>'+description[0]+'|'+description[6]+'|frameshift+2\n'+str(record.seq[2:].translate())+'\n')
    out.close()
    return 0

def frameshift_detection(input_dir,output_dir,transcriptome,template_xml):
    if not os.path.isfile(transcriptome+'-frameshift-proteome.fa'):
        print('Creating frameshift proteome from transcriptome: '+transcriptome)
        frameshift_proteome(transcriptome)
        
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
        
def noncoding_translation_detection(input_dir,output_dir,transcriptome,proteome,template_xml,noncoding_type):
    noncoding_type = noncoding_type.lower()
    path_to_custome_proteome = transcriptome+'-'+noncoding_type+'-proteome.fa'

    if not os.path.isfile(path_to_custome_proteome):
        print('Creating custom proteome from transcriptome: '+transcriptome)
        if noncoding_type == 'lncrna' or noncoding_type == 'intron':
            lncRNA_or_intron_proteome(transcriptome,noncoding_type)           
        elif noncoding_type == 'utr':
            UTR_proteome(transcriptome)
        elif noncoding_type == 'frameshift':
            frameshift_proteome(transcriptome)
        else:
            print("ERROR: unknown noncoding_type "+noncoding_type)
            exit()
    else:
        print('Using existing proteome: '+path_to_custome_proteome)
        
    if not os.path.isfile(output_dir+'/combined/txt/peptides.txt'):
        print("MaxQuant search for peptides")
        print("---------------------------------")
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        generate_xml(template_xml,input_dir,output_dir,path_to_custome_proteome)
        cmd = MaxQuantCmd+ " " + output_dir + "/mqpar.xml "
        print(cmd)
        os.system(cmd)
        os.system('wc -l '+output_dir+'/combined/txt/peptides.txt')
    else:
        print("Found existing MaxQuant search result: "+output_dir+'/combined/txt/peptides.txt')
    
    print("Removing candidate peptides also found in the normal proteome")
    npep = filter_with_proteome(output_dir+'/combined/txt/peptides.txt',proteome)
    print(str(npep)+' peptides saved to '+output_dir+'/combined/txt/peptides.txt-filtered.txt')
                  
def filter_with_proteome(path_to_pep,path_to_proteome):
    '''
    path_to_pep: path to peptides.txt
    path_to_proteome: path to reference proteome file
    '''
    try:
        pep = pd.read_csv(path_to_pep,sep = '\t',index_col='Sequence')
        allpep = open(path_to_proteome).read()
        for seq in pep.index:
            if seq in allpep:
                pep=pep.drop(seq)
        pep.to_csv(path_to_pep+'-filtered.txt')        
        return len(pep)
    except:
        return 0
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    parser.add_argument('input_dir',help='Required. Path to a folder with raw data (*.raw)')

    parser.add_argument("--output-dir",action='store',default='NA',
                         help='Output folder name. Default: same as input')

    parser.add_argument("--transcriptome", action='store',default='NA',
                         help='Path to transcriptome fasta file')
              
    parser.add_argument("--proteome", action='store',default=proteome,
                         help='Path to proteome fasta file')

    parser.add_argument("--template-xml", dest='template_xml', action='store',default=template_xml,
                         help='A template xml file')
    
    parser.add_argument("--noncoding-type", dest='noncoding_type', action='store',default='frameshift',
                         help='frameshift, UTR, lncRNA, or intron. Default is frameshift')

    '''
    parser.add_argument("--utr", action="store_true", 
                    help="run UTR peptide analysis instead of frameshift")
    
    parser.add_argument("--lncrna", action="store_true", 
                    help="run lncRNA peptide analysis instead of frameshift")
    '''
              
    args = parser.parse_args()

    # check the number of *.raw files in the input folder
    args.input_dir = os.path.abspath(args.input_dir)
    nSample = valid_raw_folder(args.input_dir)
    if nSample < 0: # folder not found
        exit()

    # the default output folder is the input folder
    if args.output_dir == 'NA':
        args.output_dir = args.input_dir
    # create the output folder if not exist
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    # create sub-folder for specific type of analysis (frameshift, utr, lncrna, intron etc)
    args.output_dir = args.output_dir +'/'+str(args.noncoding_type).lower()          
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    # default transcriptome file for each type of analysis
    if args.transcriptome == 'NA':
        if str(args.noncoding_type).lower() == 'frameshift':
            args.transcriptome  = transcriptome_frameshift
        elif str(args.noncoding_type).lower() == 'utr':
            args.transcriptome  = transcriptome_mrna
        elif str(args.noncoding_type).lower() == 'lncrna':
            args.transcriptome  = transcriptome_lncrna
        else:
            print('ERROR: unknown noncoding-type: '+args.noncoding_type)
            exit()
           
    # show parameters
    print("- noncoding type : "+args.noncoding_type)
    print("- transcriptome  : "+os.path.abspath(args.transcriptome))
    print("- proteome       : "+os.path.abspath(args.proteome))
    print("- template xml   : "+os.path.abspath(args.template_xml))
    print("- output         : "+os.path.abspath(args.output_dir))
    print("- input          : "+os.path.abspath(args.input_dir))
    print("                   " + str(nSample) + " raw file(s)")
    
    # run the analysis
    noncoding_translation_detection(args.input_dir,args.output_dir,args.transcriptome,args.proteome,args.template_xml,args.noncoding_type)
