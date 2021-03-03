'''
A pipeline for detecting frameshift events from proteomics data

Xuebing Wu

'''

# Required: please update the path to MaxQuant binary
MaxQuantCmd="dotnet /home/xw2629/software/MaxQuant/bin/MaxQuantCmd.exe"

# Option: setup default files
transcriptome_frameshift='./reference/human.CDS.fa'    # for frameshift, same as substitution
transcriptome_lncrna='./reference/human.lncRNA.fa'     # for lncRNA analysis, downloaded from GENCODE, lncRNA sequence
transcriptome_mrna='./reference/human.mRNA.fa'         # for UTR analysis, downloaded from GENCODE, protein-coding transcript sequence
transcriptome_intron='./reference/human.intron.fa'     # for intron analysis, downloaded from UCSC Table browser, gencode.v32, +9nt flanking sequence
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
        if '-' in record.seq or 'N' in record.seq:
            continue
        header=record.description.split('|') #  >ENST00000641515.2|ENSG00000186092.6|OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-202|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|
        description = '>-|'+header[0]+':'+header[5]
        txpt_len = int(header[6])
        if header[7][:3] == 'CDS':
            s,e = header[7][4:].split('-')
        elif header[8][:3] == 'CDS':
            s,e = header[8][4:].split('-')
        elif header[9][:3] == 'CDS':
            s,e = header[9][4:].split('-')
        else:
            pass
        s=int(s) # CDS start
        e=int(e) # CDS end
        # 5' UTR
        if s > 9: # at least 9nt 5' UTR
            utr5seq = record.seq[:(s+9)]
            utr5seq = utr5seq[len(utr5seq)%3:] # make it in frame with CDS by skipping the first or the first 2 bases
            #out.write(parse_fragments('>'+txpt_id+'|'+gene_symbol+'|UTR5-0|',str(utr5seq[0:].translate()),7))
            #out.write(parse_fragments('>'+txpt_id+'|'+gene_symbol+'|UTR5-1|',str(utr5seq[1:].translate()),7))
            #out.write(parse_fragments('>'+txpt_id+'|'+gene_symbol+'|UTR5-2|',str(utr5seq[2:].translate()),7))
            out.write(description+':UTR5-0|-|\n'+str(utr5seq[0:].translate())+'\n')
            out.write(description+':UTR5-1|-|\n'+str(utr5seq[1:].translate())+'\n')
            out.write(description+':UTR5-2|-|\n'+str(utr5seq[2:].translate())+'\n')
        if txpt_len - e >= 9: # 3' UTR at least 9 nt
            utr3seq = record.seq[(e-9):]
            out.write(description+':UTR3-0|-|\n'+str(utr3seq[0:].translate())+'\n')
            out.write(description+':UTR3-1|-|\n'+str(utr3seq[1:].translate())+'\n')
            out.write(description+':UTR3-2|-|\n'+str(utr3seq[2:].translate())+'\n')
    out.close()
    return 0
    
def lncRNA_or_intron_proteome(transcriptome,analysis): 
    '''
    generate protein sequences in three reading frames in lncRNAs
    - input  : a fasta file containing sequences of all lncRNAs
    - output : fasta files for proteins encoded in each reading frame, regardless of AUG
               - output in the same folder as input
               - output file suffix: -proteome.fa
    '''
    if not os.path.isfile(transcriptome):
        return -1 # file not found
    analysis = str(analysis).lower()
    
    out = open(transcriptome+'-'+analysis+'-proteome.fa','w')

    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
        record.seq = record.seq.upper()
        if '-' in record.seq or 'N' in record.seq:
            continue
        if analysis == 'lncrna': # >ENST00000326734.2|ENSG00000177757.2|OTTHUMG00000002471.2|OTTHUMT00000007025.2|FAM87B-201|FAM87B|1947|
            descriptions = record.description.split('|')
            description = '-|'+descriptions[0]+':'+descriptions[5]+':lncrna'
        else:# >hg38_wgEncodeGencodeCompV36_ENST00000641515.2_1 range=chr1:65565-69045 5'pad=9 3'pad=9 strand=+ repeatMasking=none
            descriptions = record.description.split(' ')
            description = '-|'+descriptions[0]+':'+descriptions[1]+':intron'
        #out.write(parse_fragments('>'+description+'|'+analysis+'|0|',str(record.seq[0:].translate()),7))
        #out.write(parse_fragments('>'+description+'|'+analysis+'|1|',str(record.seq[1:].translate()),7))
        #out.write(parse_fragments('>'+description+'|'+analysis+'|2|',str(record.seq[2:].translate()),7))      
        out.write('>'+description+'-0|-|\n'+str(record.seq[0:].translate())+'\n')
        out.write('>'+description+'-1|-|\n'+str(record.seq[1:].translate())+'\n')
        out.write('>'+description+'-2|-|\n'+str(record.seq[2:].translate())+'\n')
    out.close()
    return 0

# NOT USED!!!
def parse_fragments(header,translation,min_pep_len):
    '''
    generate fasta output for a translated peptides containing stop codons: 
    header: fasta header to be used
    translation: amino acid sequence containing stop codons, e.g. NNNN*NNNNN*NNNN*NN
    min_pep_len: 7, shortest peptides to consider
    
    example output:
    
    >header1
    NNNN
    >header2
    NNNNN
    >header3
    NNNN
    
    If stop codon (*) is not removed, it will match any amino acids
    
    '''
    translations = translation.split('*')
    n = 0
    fasta_block = ''
    for translation in translations:
        if len(translation) >= min_pep_len:
            n = n + 1
            fasta_block = fasta_block + header +str(n) + '\n' + translation + '\n'
    return fasta_block
    
def frameshift_proteome(transcriptome):
    '''
    generate protein sequences in three reading frames
    - input  : a fasta file containing CDS of all transcripts, each starts with start codon and ends with stop codon
    - output : a fasta file for peptides encoded in shifted reading frame. 
               - output in the same folder as input
               - output file suffix: -frameshift-proteome.fa
    '''
    # note that this doesn't include 3' UTR
    if not os.path.isfile(transcriptome):
        return -1 # file not found
    out = open(transcriptome+'-frameshift-proteome.fa','w')
    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
        record.seq = record.seq.upper()   
        if '-' in record.seq or 'N' in record.seq:
            continue
        if len(record.seq)%3 == 0 and record.seq[:3] in {'ATG','GTG','TTG','ATT','CTG'} and record.seq[-3:].translate()=='*':
            descriptions = record.description.split(' ')
            description = '>-|'+descriptions[0]+':'+descriptions[6].split(':')[1]
            #out.write(parse_fragments('>'+description[0]+'|'+description[6]+'|frameshift+1|',str(record.seq[1:].translate()),7))
            #out.write(parse_fragments('>'+description[0]+'|'+description[6]+'|frameshift+2|',str(record.seq[2:].translate()),7))
            out.write(description+':frameshift-1|-|\n'+str(record.seq[1:].translate())+'\n')
            out.write(description+':frameshift-2|-|\n'+str(record.seq[2:].translate())+'\n')
    out.close()
    return 0
        
def noncoding_translation_detection(input_dir,output_dir,transcriptome,proteome,template_xml,analysis):
    analysis = analysis.lower()
    path_to_custom_proteome = transcriptome+'-'+analysis+'-proteome.fa'

    # generate custom proteome for noncanonical analysis
    if analysis == 'canonical':
        path_to_custom_proteome = proteome
    else:
        path_to_custom_proteome = transcriptome+'-'+analysis+'-proteome.fa'
        if not os.path.isfile(path_to_custom_proteome):
            print('Creating custom proteome from transcriptome: '+transcriptome)
            if analysis == 'lncrna' or analysis == 'intron':
                lncRNA_or_intron_proteome(transcriptome,analysis)           
            elif analysis == 'utr':
                UTR_proteome(transcriptome)
            elif analysis == 'frameshift':
                frameshift_proteome(transcriptome)
            else:
                print("ERROR: unknown analysis type: "+analysis)
                exit()
        else:
            print('Using existing proteome: '+path_to_custom_proteome)
        
    path_to_evidence = output_dir+'/combined/txt/evidence.txt'
    if not os.path.isfile(path_to_evidence):
        print("MaxQuant search for peptides")
        print("---------------------------------")
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        generate_xml(template_xml,input_dir,output_dir,path_to_custom_proteome)
        cmd = MaxQuantCmd+ " " + output_dir + "/mqpar.xml "
        print(cmd)
        os.system(cmd)
        #os.system('wc -l '+output_dir+'/combined/txt/evidence.txt')
    else:
        print("Found existing MaxQuant search result: "+path_to_evidence)
    
    # filtering peptides identified 
    try:
        pep = pd.read_csv(path_to_evidence,sep = '\t',index_col='Sequence')
        print(str(len(pep))+" peptides identified")
        
        pep = pep[pep['Potential contaminant'] != '+']
        print(str(len(pep))+" peptides remain after removing potential contaminant")
        
        if analysis != 'canonical':
            allpep = open(proteome).read()
            for seq in pep.index:
                if seq in allpep:
                    pep=pep.drop(seq)
            print(str(len(pep))+" peptides remain after removing canonical peptides")
        pep.to_csv(path_to_evidence+'-filtered.txt')  
        print(str(len(pep))+' peptides saved to '+path_to_evidence+'-filtered.txt')
    except:
        print("0 peptides were identified in the MaxQuant search")
                  
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    parser.add_argument('input_dir',help='Required. Path to a folder with raw data (*.raw)')

    parser.add_argument("--analysis", dest='analysis', action='store',default='frameshift,utr,lncrna,intron',
                         help='combination of frameshift, utr, lncrna, or intron separated by comma (no space). Default: frameshift,utr,lncrna,intron')

    parser.add_argument("--output-dir",action='store',default='NA',
                         help='Output folder name. Default: same as input')

    parser.add_argument("--transcriptome", action='store',default='NA',
                         help='Path to transcriptome fasta file. See README for details')
              
    parser.add_argument("--proteome", action='store',default=proteome,
                         help='Path to proteome fasta file')

    parser.add_argument("--template-xml", dest='template_xml', action='store',default=template_xml,
                         help='A template xml file')
              
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
        
    # map to canonical proteome if not already run
    output_dir = args.output_dir +'/canonical'         
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    # show parameters
    print("analysis         : canonical ")
    print("- proteome       : "+args.proteome)
    print("- template xml   : "+args.template_xml)
    print("- output         : "+output_dir)
    print("- input          : "+args.input_dir)
    print("                   " + str(nSample) + " raw file(s)")
    noncoding_translation_detection(args.input_dir,output_dir,'',args.proteome,args.template_xml,'canonical')
        
    # custom proteome analysis
    analyses = args.analysis.split(',')
    custom_transcriptome = ''
    for analysis in analyses:
        analysis = analysis.lower()
        # default transcriptome file for each type of analysis
        if args.transcriptome == 'NA':
            if analysis == 'frameshift':
                custom_transcriptome  = transcriptome_frameshift
            elif analysis == 'utr':
                custom_transcriptome  = transcriptome_mrna
            elif analysis == 'lncrna':
                custom_transcriptome  = transcriptome_lncrna
            elif analysis == 'intron':
                custom_transcriptome  = transcriptome_intron
            else:
                print('ERROR: unknown noncoding-type: '+analysis)
                exit()
        
        # create sub-folder for specific type of analysis (frameshift, utr, lncrna, intron etc)
        output_dir = args.output_dir +'/'+analysis         
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        
        # show parameters
        print("analysis         : "+analysis)
        print("- transcriptome  : "+custom_transcriptome)
        print("- proteome       : "+args.proteome)
        print("- template xml   : "+args.template_xml)
        print("- output         : "+output_dir)
        print("- input          : "+args.input_dir)
        print("                   " + str(nSample) + " raw file(s)")

        # run the analysis
        noncoding_translation_detection(args.input_dir,output_dir,custom_transcriptome,args.proteome,args.template_xml,analysis)

