'''
A pipeline for detecting translation errors from proteomics data

https://github.com/xuebingwu/tfams

Xuebing Wu
'''

# Required: please update the path to MaxQuant binary
MaxQuantCmd="dotnet /home/xw2629/software/MaxQuant/bin/MaxQuantCmd.exe"

### Optional but required if want to use default settings

# default of --standard-xml
template_xml_standard='./template_xml/mqpar-standard.xml'         # Standard MaxQuant search parameters, no dependent peptide search

# default of --substitution-xml
template_xml_substitution='./template_xml/mqpar-substitution.xml' # MaxQuant parameters with dependent peptide search and match between runs

# default of --proteome
proteome='./reference/human.protein.fa'                # amino acid sequences of all proteins

# default of --transcriptome
transcriptome='./reference/human.CDS.fa'               # for substitution, only CDS of mRNAs

# default for other transcriptomes used to generate noncanonical peptide databases 
transcriptome_frameshift='./reference/human.CDS.fa'    # for frameshift analysis, same as substitution
transcriptome_lncrna='./reference/human.lncRNA.fa'     # for lncRNA analysis, downloaded from GENCODE, lncRNA sequence
transcriptome_mrna='./reference/human.mRNA.fa'         # for UTR analysis, downloaded from GENCODE, protein-coding transcript sequence
transcriptome_intron='./reference/human.intron.fa'     # for intron analysis, downloaded from UCSC Table browser, gencode.v32, +9nt flanking sequence

import argparse
import os
from generate_xml import *
from Bio import SeqIO
import time
import pandas as pd

def generate_UTR_proteome(transcriptome):
    '''
    input: 
        - a fasta file of mRNA sequences, including UTRs. 
        - downloaded from GENCODE, sequenes for protein-coding transcripts
        - header format: >ENST00000641515.2|ENSG00000186092.6|OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-202|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|\
    output:
        - a fasta file of amino acid sequences resulted from translating 5' UTRs or 3' UTRs 
        - output file name: *-utr-proteome.fa, in the same folder as input
        - allow translation in all 3 reading frames
        - allow stop codons to remain in the peptide sequences (as *)
        - header format: >en|ENST00000641515.2_OR4F5_UTR5_0|
        - UTR5_0 means translation of 5' UTR in frame 0, same as the coding sequence downstream
        - UTR3_0 means translation of 3' UTR in frame 0, same as the coding sequence upstream
    '''
    if not os.path.isfile(transcriptome):
        return -1 # file not found
    out = open(transcriptome+'-utr-proteome.fa','w')
    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
        record.seq = record.seq.upper()
        if '-' in record.seq or 'N' in record.seq:
            continue
        header=record.description.split('|') 
        description = '>en|'+header[0]+'_'+header[5]
        txpt_len = int(header[6])
        if header[7][:3] == 'CDS':
            s,e = header[7][4:].split('-')
        elif header[8][:3] == 'CDS':
            s,e = header[8][4:].split('-')
        elif header[9][:3] == 'CDS':
            s,e = header[9][4:].split('-')
        else:
            continue
        s=int(s) # CDS start
        e=int(e) # CDS end
        # 5' UTR
        if s > 9: # at least 9nt 5' UTR
            utr5seq = record.seq[:(s+9)]
            utr5seq = utr5seq[len(utr5seq)%3:] # make it in frame with CDS by skipping the first or the first 2 bases
            out.write(description+'_UTR5_0|\n'+str(utr5seq[0:].translate())+'\n')
            out.write(description+'_UTR5-1|\n'+str(utr5seq[1:].translate())+'\n')
            out.write(description+'_UTR5-2|\n'+str(utr5seq[2:].translate())+'\n')
        if txpt_len - e >= 9: # 3' UTR at least 9 nt
            utr3seq = record.seq[(e-9):]
            out.write(description+'_UTR3_0|\n'+str(utr3seq[0:].translate())+'\n')
            out.write(description+'_UTR3_1|\n'+str(utr3seq[1:].translate())+'\n')
            out.write(description+'_UTR3_2|\n'+str(utr3seq[2:].translate())+'\n')
    out.close()
    return 0
    
def generate_lncRNA_or_intron_proteome(transcriptome,analysis): 
    '''
    input
        - a fasta file containing sequences of all lncRNAs or introns
        - lncRNA: downloaded from GENCODE, sequences of all lncRNAs
        - lncRNA header: >ENST00000326734.2|ENSG00000177757.2|OTTHUMG00000002471.2|OTTHUMT00000007025.2|FAM87B-201|FAM87B|1947|
        - intron: downloaded from UCSC Table browser, with 9nt flanking sequence, by chromosome and then merged together
    output
        - a fasta files for proteins encoded in each reading frame, regardless of start or stop codons
        - output in the same folder as input, file name: *-lncrna-proteome.fa or *-intron-proteome.fa
        - header for lncRNA: >en|ENST00000326734.2_FAM87B_lncrna_0 
        - header for itnron: >en|ENST00000326734.2_range=chr1:65565-69045_intron_0
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
            description = '>en|'+descriptions[0]+'_'+descriptions[5]+'_lncrna_'
        else:# >hg38_wgEncodeGencodeCompV36_ENST00000641515.2_1 range=chr1:65565-69045 5'pad=9 3'pad=9 strand=+ repeatMasking=none
            descriptions = record.description.split(' ')
            description = '>en|'+descriptions[0].split('_')[2]+'_'+descriptions[1]+'_intron_'    
        out.write(description+'0|\n'+str(record.seq[0:].translate())+'\n')
        out.write(description+'1|\n'+str(record.seq[1:].translate())+'\n')
        out.write(description+'2|\n'+str(record.seq[2:].translate())+'\n')
    out.close()
    return 0
    
def generate_frameshift_proteome(transcriptome):
    '''
    generate protein sequences in three reading frames
    input: 
        - a fasta file containing CDS of all transcripts, each starts with start codon and ends with stop codon
        - noncanonical start codon allowed ('GTG','TTG','ATT','CTG')
        - downloaded from emsembl http://ftp.ensembl.org/pub/release-103/fasta/
        - header: >ENST00000632684.1 cds chromosome:GRCh38:7:142786213:142786224:1 gene:ENSG00000282431.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:TRBD1 description:T cell receptor beta diversity 1 [Source:HGNC Symbol;Acc:HGNC:12158]
    output:
        - a fasta file for peptides encoded in shifted reading frame (1 or 2). 
        - output in the same folder as input, file name: *-frameshift-proteome.fa
        - header: >en|ENST00000632684.1_TRBD1_frameshift_1
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
            gene_symbol = descriptions[0]+'_'+descriptions[6].split(':')[1]+'_frameshift_1'
            description = 'en|'+gene_symbol+'|'
            out.write('>'+description+'\n'+str(record.seq[1:].translate())+'\n')
            gene_symbol = descriptions[0]+'_'+descriptions[6].split(':')[1]+'_frameshift_2'
            description = 'en|'+gene_symbol+'|'
            out.write('>'+description+'\n'+str(record.seq[2:].translate())+'\n')
    out.close()
    return 0
        
def generate_custom_proteome(analysis):
    analysis = analysis.lower()
    
    print('Creating custom proteome from transcriptome: ')
    if analysis == 'lncrna':
        path_to_custom_proteome = transcriptome_lncrna+'-lncrna-proteome.fa'
        if not os.path.isfile(path_to_custom_proteome):
            generate_lncRNA_or_intron_proteome(transcriptome_lncrna,'lncrna') 
        else:
            print('Custom proteome already exists. Please delete it to create a new one: '+path_to_custom_proteome)
    elif analysis == 'intron':
        path_to_custom_proteome = transcriptome_intron+'-intron-proteome.fa'
        if not os.path.isfile(path_to_custom_proteome):
            generate_lncRNA_or_intron_proteome(transcriptome_intron,'intron') 
        else:
            print('Custom proteome already exists. Please delete it to create a new one: '+path_to_custom_proteome)
    elif analysis == 'utr':
        path_to_custom_proteome = transcriptome_mrna+'-utr-proteome.fa'
        if not os.path.isfile(path_to_custom_proteome):
            generate_UTR_proteome(transcriptome_mrna)
        else:
            print('Custom proteome already exists. Please delete it to create a new one: '+path_to_custom_proteome)
    elif analysis == 'frameshift':
        path_to_custom_proteome = transcriptome_frameshift+'-frameshift-proteome.fa'
        if not os.path.isfile(path_to_custom_proteome):
            generate_frameshift_proteome(transcriptome_frameshift)
        else:
            print('Custom proteome already exists. Please delete it to create a new one: '+path_to_custom_proteome)
    else:
        print("ERROR: unknown analysis type: "+analysis)
        exit()
        
    return path_to_custom_proteome

def maxquant_standard_search(path_to_proteome,input_dir,output_dir,template_xml):
    # standard maxquant analysis
    
    # check if already done
    path_to_evidence = output_dir+'/combined/txt/evidence.txt'
    if not os.path.isfile(path_to_evidence):
        print("MaxQuant search for peptides")
        print("---------------------------------")
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        generate_xml(template_xml,input_dir,output_dir,path_to_proteome)
        cmd = MaxQuantCmd+ " " + output_dir + "/mqpar.xml "
        print(cmd)
        os.system(cmd)
    else:
        print("MaxQuant result found in the output folder. To run again please delete existing output folder: "+output_dir)
    
def filter_maxquant_result(output_dir,path_to_proteome):
    '''
    filter 1: potential contamiants
    filter 2: peptides found in the normal/canonical proteome
    
    input:
        - output_dir: a parental folder where 'combined' from a MaxQuant run can be found
        - path_to_proteome: a fasta file with all canonical protein sequences
    output:
        - shown on screen how many peptides remain after each filter
        - filtered peptides saved to output_dir+'/combined/txt/evidence.txt-filtered.txt'
    '''
    path_to_evidence = output_dir+'/combined/txt/evidence.txt'
    try:
        pep = pd.read_csv(path_to_evidence,sep = '\t')
        print(str(len(pep))+" peptides identified")
    except:
        print("0 peptides were identified in the MaxQuant search")
        
    pep = pep[pep['Potential contaminant'] != '+']
    print(str(len(pep))+" peptides remain after removing potential contaminant")
        
    # remove peptides also found in canonical proteome
    # note that two sequences can be the same, so do not index using sequence
    allpep = open(path_to_proteome).read()
    for i in pep.index:
        #print(pep.at[i,'Sequence'])
        if pep.at[i,'Sequence'] in allpep:
            pep=pep.drop(i)
    print(str(len(pep))+" peptides remain after removing canonical peptides")
    pep.to_csv(path_to_evidence+'-filtered.txt')  
    print('filtered peptides saved to '+path_to_evidence+'-filtered.txt')        
    
# NOT USED!!! slows the analysis by a lot
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
    '''
    translations = translation.split('*')
    n = 0
    fasta_block = ''
    for translation in translations:
        if len(translation) >= min_pep_len:
            n = n + 1
            fasta_block = fasta_block + header +str(n) + '\n' + translation + '\n'
    return fasta_block

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    parser.add_argument('input_dir',help='Required. Path to a folder with raw data (*.raw)')

    parser.add_argument("--analysis", dest='analysis', action='store',default='NA',
                         help='combination of frameshift, utr, lncrna, intron, and substitution separated by comma (no space)')

    parser.add_argument("--output-dir",action='store',default='NA',
                         help='Output folder name. Default: same as input')

    parser.add_argument("--transcriptome", action='store',default=transcriptome,
                         help='Path to transcriptome fasta file. See README for details')
              
    parser.add_argument("--proteome", action='store',default=proteome,
                         help='Path to proteome fasta file')

    parser.add_argument("--substitution-xml", dest='template_xml_substitution', action='store',default=template_xml_substitution,
                         help='A template xml file for substitution detection')
    
    parser.add_argument("--standard-xml", dest='template_xml_standard', action='store',default=template_xml_standard,
                         help='A template xml file for standard MaxQuant search')
              
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
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    # show parameters
    print("analysis         : standard ")
    print("- proteome       : "+args.proteome)
    print("- template xml   : "+args.template_xml_standard)
    print("- output         : "+args.output_dir)
    print("- input          : "+args.input_dir)
    print("                   " + str(nSample) + " raw file(s)")
    maxquant_standard_search(args.proteome,args.input_dir,args.output_dir,args.template_xml_standard)

    # custom proteome analysis, if specified by --analysis
    if args.analysis == 'NA':
        exit()

    analyses = args.analysis.split(',')

    for analysis in analyses:

        analysis = analysis.lower()
        
        if analysis == 'substitution':
            continue
        
        path_to_custom_proteome = generate_custom_proteome(analysis)

        # create sub-folder for specific type of analysis (frameshift, utr, lncrna, intron etc)
        output_dir = args.output_dir +'/'+analysis         
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        # show parameters
        print("analysis         : "+analysis)
        print("- proteome       : "+path_to_custom_proteome)
        print("- template xml   : "+args.template_xml_standard)
        print("- output         : "+output_dir)
        print("- input          : "+args.input_dir)
        print("                   " + str(nSample) + " raw file(s)")

        # run the analysis
        maxquant_standard_search(path_to_custom_proteome,args.input_dir,output_dir,args.template_xml_standard)
        filter_maxquant_result(output_dir,args.proteome)
        
    if 'substitution' in args.analysis:
        
        args.output_dir = args.output_dir +'/substitution'
        
        # create params.py
        f = open("params.py", "w")
        f.write("input_dir = '"+args.input_dir+"'\n")
        f.write("output_dir = '"+args.output_dir+"'\n")
        f.write("transcriptome = '"+args.transcriptome+"'\n")
        f.close()    
        os.system('cat params.core >> params.py')

        print("Path to files:")
        print("- proteome          : "+args.proteome)
        print("- transcriptome     : "+args.transcriptome)
        print("- template xml      : "+args.template_xml_substitution)
        print("- output            : "+args.output_dir)
        print("- input             : "+args.input_dir)
        print("                      " + str(nSample) + " raw file(s)")

        if os.path.isfile(args.output_dir+'/combined/txt/allPeptides.txt'):
            print("MaxQuant search results detected in the output folder "+args.output_dir)
            print("- to run MaxQuant again, please change output directory")
        else:
            print("Generate MaxQuant parameter file (xml)")
            generate_xml(args.template_xml_substitution,args.input_dir,args.output_dir,args.proteome)
            print("MaxQuant search of dependent peptides")
            cmd = MaxQuantCmd+ " " + args.output_dir + "/mqpar.xml"
            print('- '+cmd)
            os.system(cmd)
        print("Detection and filtering")
        import detect
        import quantify
        import plot