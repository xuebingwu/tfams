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

# default of --variant
variant='./reference/human.variant.fa'                 # peptides encoded by SNV

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
    out0 = open(transcriptome+'-proteome.fa','w')
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
            # generate non-frameshifted proteins
            gene_symbol = descriptions[0]+'_'+descriptions[6].split(':')[1]+'_frameshift_0'
            description = 'en|'+gene_symbol+'|'
            out0.write('>'+description+'\n'+str(record.seq.translate())+'\n')
    out.close()
    out0.close()
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
        print("Deleting intermediate files from MaxQuant run")
        clean_up(input_dir,output_dir)
    else:
        print("MaxQuant result found in the output folder. To run again please delete existing output folder: "+output_dir)
    
def filter_maxquant_result(output_dir,path_to_proteome,path_to_variant_peptide):
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
    pep = pep[pep['Reverse'] != '+']
    print(str(len(pep))+" peptides remain after removing potential contaminant and reverse match")
        
    # remove peptides also found in canonical proteome
    # note that two sequences can be the same, so do not index using sequence
    allpep = open(path_to_proteome).read().replace("\n","")
    for i in pep.index:
        #print(pep.at[i,'Sequence'])
        if pep.at[i,'Sequence'] in allpep:
            pep=pep.drop(i)
    print(str(len(pep))+" peptides remain after removing canonical peptides")
    
    pep = pep[~pep['Intensity'].isna()]
    pep = pep[pep['Intensity']>0]
    print(str(len(pep))+" peptides remain after removing peptides with no positive intensity values")

    pep = pep[~pep['Proteins'].isna()]
    print(str(len(pep))+" peptides remain after removing peptides with no assigned proteins")
    
    # SNPs
    pep['SNV'] = 0
    n=0
    if os.path.isfile(path_to_variant_peptide):
        variant_pep = open(path_to_variant_peptide).read().replace("\n","")
        for i in pep.index:
            if pep.at[i,'Sequence'] in variant_pep:
                #pep=pep.drop(i)
                pep.at[i,'SNV'] = 1
                n=n+1
        print(str(n)+" peptides are marked as potential SNP peptides")
    else:
        print("Skip variant filter: variant peptide file not found or not set: "+path_to_variant_peptide)
    
    pep.to_csv(path_to_evidence+'-filtered.txt')  
    print('filtered peptides saved to '+path_to_evidence+'-filtered.txt')
    return len(pep)
    
def quantification_and_plot(output_dir,canonical_output_dir):
    # count, median, mean intensity of identified peptides
    # count, median, mean intensity of protein groups with identified peptides
    # same for canonical output
    # violin plot for the distribution of peptide/protein intensities
    
    #output_dir = '/home/xw2629/xuebing/amino-acid-substitution/raw-data/frameshift-test/frameshift/'
    #output_dir = '/home/xw2629/xuebing/amino-acid-substitution/raw-data/frameshift-test/canonical/'

    pep = pd.read_csv(output_dir+'/combined/txt/evidence.txt-filtered.txt')
    
    pep = pep[~pep['Intensity'].isna() & pep['Intensity']>0]
    
    print(" peptides identified: "+str(len(pep)))
    print(" peptide intensity, median/1e6: "+str(pep['Intensity'].median()/1e6))
    print(" peptide intensity, mean/1e6: "+str(pep['Intensity'].mean()/1e6))
    #plt.hist(np.log10(pep['Intensity']),bins=100)

    plt.violinplot([np.log10(pep['Intensity'])],showmeans=False,showmedians=True)
    plt.ylabel('log10 (Intensity)')
    plt.xticks([1], ['Pep'])
    plt.show()
    
    # protein level
    protein_intensity = {}
    for i in pep.index:
        protein = pep.loc[i]['Leading razor protein'].split('_')[1] # ENST00000564104.5_TPPP3_frameshift_2
        if not (protein in protein_intensity):
            protein_intensity[protein] = 0
        protein_intensity[protein] = protein_intensity[protein] + pep.loc[i]['Intensity']
    protein_intensity_df = pd.DataFrame.from_dict(protein_intensity, orient='index',columns=['Intensity'])
    print(" proteins identified: "+str(len(protein_intensity_df)))
    print(" protein intensity, median/1e6: "+str(protein_intensity_df.Intensity.median()/1e6))
    print(" peptide intensity, mean/1e6: "+str(protein_intensity_df.Intensity.mean()/1e6))
    
    plt.violinplot([np.log10(protein_intensity_df.Intensity)],showmeans=False,showmedians=True)
    plt.ylabel('log10 (Intensity)')
    plt.xticks([1], ['Protein'])

def clean_up(input_dir,output_dir):
    '''
    input_dir: for each X.raw file, remove corresponding X.index and folder X
    output_dir: 
    '''
    
    # make sure output file exists
    if not os.path.isfile(output_dir+'/combined/txt/evidence.txt'):
        print("warning: it appears the analysis has not been completed. skip clean up")
        return 0
    
    # delete everything under 'combined' except 'txt'
    os.system('mv '+output_dir+'/combined/txt '+output_dir+'/')
    os.system('rm -rf '+output_dir+'/combined/*')
    os.system('mv '+output_dir+'/txt '+output_dir+'/combined/')
    
    # clean up intermediate files in input folder 
    if os.path.isdir(input_dir):
        files = os.listdir(input_dir)
        for file in files:
            if file[-4:] == ".raw":
                os.system('rm '+input_dir+'/'+file[:-4]+'.index')
                os.system('rm -rf '+input_dir+'/'+file[:-4])
    else:
        print("warning: input folder not found: "+input_dir)
    
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

    parser.add_argument("--analysis", dest='analysis', action='store',default='canonical',
                         help='combination of canonical, frameshift, utr, lncrna, intron, and substitution separated by comma (no space). Default: canonical')

    parser.add_argument("--output-dir",action='store',default='NA',
                         help='Output folder name. Default: same as input')

    parser.add_argument("--transcriptome", action='store',default=transcriptome,
                         help='Path to transcriptome fasta file. See README for details')
              
    parser.add_argument("--proteome", action='store',default=proteome,
                         help='Path to proteome fasta file')
    
    parser.add_argument("--variant", action='store',default=variant,
                         help='Path to variant peptides caused by SNV')

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

    analyses = args.analysis.split(',')

    for analysis in analyses:

        analysis = analysis.lower()
        
        if analysis == 'substitution':
            continue
        
        path_to_custom_proteome = args.proteome
        if analysis != 'canonical':
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
        
        if analysis != 'canonical':
            npep = filter_maxquant_result(output_dir,args.proteome,args.variant)
        
    if 'substitution' in args.analysis:
        
        args.output_dir = args.output_dir +'/substitution'
        
        # create params.py
        f = open("params.py", "w")
        f.write("input_dir = '"+args.input_dir+"'\n")
        f.write("output_dir = '"+args.output_dir+"'\n")
        f.write("transcriptome = '"+args.transcriptome+"'\n")
        f.write("path_to_variant_peptide = '"+args.variant+"'\n")
        f.close()    
        os.system('cat params.core >> params.py')

        print("Path to files:")
        print("- proteome          : "+args.proteome)
        print("- variant           : "+args.variant)
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
            print("Deleting intermediate files from MaxQuant run")
            clean_up(args.input_dir,args.output_dir)
        print("Detection and filtering")
        import detect
        import quantify
        import plot