'''
A pipeline for detecting translation errors from proteomics data

https://github.com/xuebingwu/tfams

Xuebing Wu

'''

# Required: please update the path to MaxQuant binary
MaxQuantCmd="dotnet /home/xw2629/software/MaxQuant/bin/MaxQuantCmd.exe"

### Optional but required if default settings will be used

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

# default of --contaminant
contaminant='./reference/contaminant.fa'               # for substitution, from MaxQuant

# default for other transcriptomes used to generate noncanonical peptide databases 
transcriptome_frameshift='./reference/human.CDS.fa'    # for frameshift analysis, same as substitution
transcriptome_lncrna='./reference/human.lncRNA.fa'     # for lncRNA analysis, downloaded from GENCODE, lncRNA sequence
transcriptome_mrna='./reference/human.mRNA.fa'         # for UTR analysis, downloaded from GENCODE, protein-coding transcript sequence
transcriptome_intron='./reference/human.intron.fa'     # for intron analysis, downloaded from UCSC Table browser, gencode.v32, +9nt flanking sequence

import argparse
import os
from generate_xml import *
from substitution import *
from Bio import SeqIO
import time
import pandas as pd
import numpy as np
import random


# for plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# noncanonical translation

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
    
def generate_lncRNA_proteome(transcriptome,analysis): 
    '''
    input
        - a fasta file containing sequences of all lncRNAs
        - lncRNA: downloaded from GENCODE, sequences of all lncRNAs
        - lncRNA header: >ENST00000326734.2|ENSG00000177757.2|OTTHUMG00000002471.2|OTTHUMT00000007025.2|FAM87B-201|FAM87B|1947|
    output
        - a fasta files for proteins encoded in each reading frame, regardless of start or stop codons
        - output in the same folder as input, file name: *-lncrna-proteome.fa or *-intron-proteome.fa
        - header for lncRNA: >en|ENST00000326734.2_FAM87B_lncrna_0 
    '''
    if not os.path.isfile(transcriptome):
        return -1 # file not found
    
    out = open(transcriptome+'-'+analysis+'-proteome.fa','w')

    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
        record.seq = record.seq.upper()
        if '-' in record.seq or 'N' in record.seq:
            continue
        descriptions = record.description.split('|')
        description = '>en|'+descriptions[0]+'_'+descriptions[5]+'_lncrna_'
        out.write(description+'0|\n'+str(record.seq[0:].translate())+'\n')
        out.write(description+'1|\n'+str(record.seq[1:].translate())+'\n')
        out.write(description+'2|\n'+str(record.seq[2:].translate())+'\n')
    out.close()
    return 0

def generate_intron_proteome(transcriptome,analysis): 
    '''
    given introns are huge, stop at the first stop codon for each frame
    
    input
        - a fasta file containing sequences of all introns
        - intron: downloaded from UCSC Table browser, with 9nt flanking sequence, by chromosome and then merged together
    output
        - a fasta files for proteins encoded in each reading frame, regardless of start or stop codons
        - output in the same folder as input, file name: *-intron-proteome.fa
        - header for itnron: >en|ENST00000326734.2_range=chr1:65565-69045_intron_0
    '''
    if not os.path.isfile(transcriptome):
        return -1 # file not found
    
    out = open(transcriptome+'-'+analysis+'-proteome.fa','w')

    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
        record.seq = record.seq.upper()
        if '-' in record.seq or 'N' in record.seq:
            continue
        # >hg38_wgEncodeGencodeCompV36_ENST00000641515.2_1 range=chr1:65565-69045 5'pad=9 3'pad=9 strand=+ repeatMasking=none
        descriptions = record.description.split(' ')
        description = '>en|'+descriptions[0].split('_')[2]+'_'+descriptions[1]+'_intron_' 
        seq = str(record.seq[0:].translate().split('*')[0])
        if len(seq) > 6:
            out.write(description+'0|\n'+seq+'\n')
        seq = str(record.seq[1:].translate().split('*')[0])
        if len(seq) > 6:
            out.write(description+'1|\n'+seq+'\n')
        seq = str(record.seq[2:].translate().split('*')[0])
        if len(seq) > 6:
            out.write(description+'2|\n'+seq+'\n')
    out.close()
    return 0
    
def shuffle(seq):
    str_var = list(seq)
    random.shuffle(str_var)
    return ''.join(str_var)

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
    out1 = open(transcriptome+'-frameshift-shuffle-proteome.fa','w')
    for record in SeqIO.parse(open(transcriptome,'r'),'fasta'):
        record.seq = record.seq.upper()   
        if '-' in record.seq or 'N' in record.seq:
            continue
        if len(record.seq)%3 == 0 and record.seq[:3] in {'ATG','GTG','TTG','ATT','CTG'} and record.seq[-3:].translate()=='*':
            descriptions = record.description.split(' ')
            gene_symbol = descriptions[0]+'_'+descriptions[6].split(':')[1]+'_frameshift_1'
            description = 'en|'+gene_symbol+'|'
            seq = str(record.seq[1:].translate())
            out.write('>'+description+'\n'+seq+'\n')
            
            # shuffled frame 1
            out1.write('>'+description+'\n'+shuffle(seq)+'\n')
            
            gene_symbol = descriptions[0]+'_'+descriptions[6].split(':')[1]+'_frameshift_2'
            description = 'en|'+gene_symbol+'|'
            
            seq = str(record.seq[2:].translate())
            out.write('>'+description+'\n'+seq+'\n')
            # shuffled frame 2
            out1.write('>'+description+'\n'+shuffle(seq)+'\n')
            
            # generate non-frameshifted proteins
            gene_symbol = descriptions[0]+'_'+descriptions[6].split(':')[1]+'_frameshift_0'
            description = 'en|'+gene_symbol+'|'
            out0.write('>'+description+'\n'+str(record.seq.translate())+'\n')
    out.close()
    out0.close()
    out1.close()
    return 0
        
def generate_custom_proteome(analysis):
    
    print('Creating custom proteome from transcriptome: ')
    if analysis == 'lncrna':
        path_to_custom_proteome = transcriptome_lncrna+'-lncrna-proteome.fa'
        if not os.path.isfile(path_to_custom_proteome):
            generate_lncRNA_proteome(transcriptome_lncrna,'lncrna') 
        else:
            print('Custom proteome already exists. Please delete it to create a new one: '+path_to_custom_proteome)
    elif analysis == 'intron':
        path_to_custom_proteome = transcriptome_intron+'-intron-proteome.fa'
        if not os.path.isfile(path_to_custom_proteome):
            generate_intron_proteome(transcriptome_intron,'intron') 
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


# maxquant 

def maxquant_standard_search(path_to_proteome,input_dir,output_dir,template_xml,thread):
    # standard maxquant analysis
    
    # check if already done
    path_to_evidence = output_dir+'/combined/txt/evidence.txt'
    if not os.path.isfile(path_to_evidence):
        print("MaxQuant search for peptides")
        print("---------------------------------")
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        generate_xml(template_xml,input_dir,output_dir,path_to_proteome,thread)
        cmd = MaxQuantCmd+ " " + output_dir + "/mqpar.xml "
        print(cmd)
        os.system(cmd)
        print("Deleting intermediate files from MaxQuant run")
        clean_up_after_maxquant_run(input_dir,output_dir)
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
    print(str(len(pep))+" peptides remain after removing potential contaminant")
    pep = pep[pep['Reverse'] != '+']
    print(str(len(pep))+" peptides remain after removing reverse match")
        
    # remove peptides also found in canonical proteome
    # note that two sequences can be the same, so do not index using sequence
    allpep = open(path_to_proteome).read().replace("\n","")
    for i in pep.index:
        if pep is None:
            print("0 peptides remain after removing canonical peptides")
            return 0
            break
        if pep.at[i,'Sequence'] in allpep:
            pep=pep.drop(i,inplace=True)
    
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
    
'''
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
'''


def clean_up_after_maxquant_run(input_dir,output_dir):
    '''
    input_dir: for each X.raw file, remove corresponding X.index and folder X
    output_dir: 
    '''
    
    # make sure output file exists
    if not os.path.isfile(output_dir+'/combined/txt/allPeptides.txt'):
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
    
def substitution_second_pass(path_to_subs,input_dir,output_dir,standard_xml,thread):
    '''
    create a peptide database with only substitutions and their base peptide
    then run maxquant standard search again to get quantifications
    '''
    
    path_to_2nd_evidence = output_dir+'/second_pass/combined/txt/evidence.txt'
    
    if not os.path.isfile(path_to_2nd_evidence):
        
        subs = pd.read_csv(path_to_subs)

        print("- 2nd pass: get unique sequences for each sub and base peptide")
        bp_seqs = {}
        for i in subs.index:
            bp_seqs[subs.at[i,'modified_sequence']] = subs.at[i,'DP Base Sequence']
        print('- 2nd pass: '+str(len(bp_seqs))+' unique peptides with substitutions')
        print("- 2nd pass: write into fasta file: "+path_to_subs+'-uniq-seq.fa')
        out = open(path_to_subs+'-uniq-seq.fa','w')
        for key in bp_seqs:
            out.write('>sub|'+key+'|\n'+key+'\n')
            out.write('>sub|'+bp_seqs[key]+'|\n'+bp_seqs[key]+'\n')
        out.close()

        print("- 2nd pass: maxquant search")
        maxquant_standard_search(path_to_subs+'-uniq-seq.fa',input_dir,output_dir+'/second_pass',standard_xml,thread)
    else:
        print("- 2nd pass already done. To run again, please delete "+path_to_2nd_evidence)
    
def substitution_quantification_with_2nd_pass(path_to_subs,path_to_2nd_evidence,output_dir):
    subs = pd.read_csv(path_to_subs)
    evidence2 = pd.read_csv(path_to_2nd_evidence,sep='\t')
    intensities = {}
    min_PEP_intensity = {} # min PEP
    min_PEP = {} # min PEP
    
    # add up the intensity values for the same peptides
    for i in evidence2.index:
        if evidence2.at[i,'Intensity'] > 0:
            seq = evidence2.at[i,'Sequence']
            if not (seq in intensities):
                intensities[seq] = 0
                min_PEP_intensity[seq] = 0
                min_PEP[seq] = 1
            intensities[seq] += evidence2.at[i,'Intensity']
            if evidence2.at[i,'PEP'] < min_PEP[seq]:
                min_PEP[seq] = evidence2.at[i,'PEP']
                min_PEP_intensity[seq] = evidence2.at[i,'Intensity']
    
    # calculate intensity ratio
    uniq_subs = pd.DataFrame([],columns=['modified_sequence',
                                         'reference_sequence',
                                         'modified_intensity',
                                         'reference_intensity',
                                         'log10_mod_to_ref_intensity_ratio',
                                         'modified_min_PEP',
                                         'reference_min_PEP',
                                         'modified_intensity_min_PEP',
                                         'reference_intensity_min_PEP',
                                         'log10_mod_to_ref_intensity_ratio_min_PEP',
                                         'longest_protein',
                                         'protein_length',
                                         'position',
                                         'position_relative',
                                         'longest_protein_status',
                                         'codon',
                                         'destination',
                                         'origin',
                                         'mispairing',
                                         'positions',
                                         'proteins',
                                         'coding_status',
                                         'Retention time',
                                         'DP Time Difference'])
            
    computed=[]
    for i in subs.index:
        dpseq = subs.at[i,'modified_sequence']
        bpseq = subs.at[i,'DP Base Sequence']
        if not(dpseq in computed):
            computed.append(dpseq)
            if dpseq in intensities and bpseq in intensities:
                uniq_subs=uniq_subs.append(pd.DataFrame([[dpseq,
                                                          bpseq,
                                                          intensities[dpseq],
                                                          intensities[bpseq],
                                                          np.log10(intensities[dpseq]/intensities[bpseq]),
                                                          min_PEP[dpseq],
                                                          min_PEP[bpseq],
                                                          min_PEP_intensity[dpseq],
                                                          min_PEP_intensity[bpseq],
                                                          np.log10(min_PEP_intensity[dpseq] / min_PEP_intensity[bpseq]),
                                                          subs.at[i,'longest_protein'],
                                                          subs.at[i,'protein_length'],
                                                          subs.at[i,'position'],
                                                          subs.at[i,'position_relative'],
                                                          subs.at[i,'longest_protein_status'],
                                                          subs.at[i,'codon'],
                                                          subs.at[i,'destination'],
                                                          subs.at[i,'origin'],
                                                          subs.at[i,'mispairing'],
                                                          subs.at[i,'positions'],
                                                          subs.at[i,'proteins'],
                                                          subs.at[i,'coding_status'],
                                                          subs.at[i,'Retention time'],
                                                          subs.at[i,'DP Time Difference']
                                                         ]], columns=uniq_subs.columns))
    uniq_subs.to_csv(os.path.join(output_dir,'uniq_subs.csv'))
    
def substitution_filter(path_to_uniq_subs,path_to_variant_peptide,path_to_proteome,path_to_contaminant):
    
    '''
    if os.path.isfile(path_to_uniq_subs[:-4]+'-filtered.csv'):
        print("- skip variant filter: results already found at: "+path_to_uniq_subs[:-4]+'-filtered.csv')
        return 0
    '''
        
    uniq_subs = pd.read_csv(path_to_uniq_subs)
    print("- unique peptides identified: "+str(len(uniq_subs)))
    
    # remove peptides whose base peptides cannot be mapped to the proteome, mostly trypsin peptides
    uniq_subs = uniq_subs[pd.notnull(uniq_subs['proteins'])]
    print("- after removing those with base peptides cannot be mapped to the proteome : "+str(len(uniq_subs)))
    
    # remove peptides also found in canonical proteome
    n=0
    if os.path.isfile(path_to_proteome):
        cont = open(path_to_proteome).read().replace("\n","")
        for i in uniq_subs.index:
            if uniq_subs.at[i,'modified_sequence'] in cont:
                uniq_subs.drop(i,inplace=True)
                n=n+1
        print("- peptides discarded for matching canonical proteome: "+str(n))
    else:
        print("- skip proteome filter: file not found or not set: "+path_to_proteome)
        
    # remove peptides also found in potential contaminant sequence
    #uniq_subs['contaminant'] = 0
    n=0
    if os.path.isfile(path_to_contaminant):
        cont = open(path_to_contaminant).read().replace("\n","")
        for i in uniq_subs.index:
            if uniq_subs.at[i,'modified_sequence'] in cont:
                #uniq_subs.at[i,'contaminant'] = 1
                uniq_subs.drop(i,inplace=True)
                n=n+1
        print("- peptides discarded as potential contaminant: "+str(n))
    else:
        print("- skip contaminant filter: file not found or not set: "+path_to_contaminant)

    print("- likely antibody sequences: "+str(sum(uniq_subs['longest_protein_status'] != 'normal')))
    
    # mark substitutions that have been observed before, in recombinant proteins (usually CHO cells) or native proteins in mammalian cells
    # https://www.sciencedirect.com/science/article/pii/S073497501730126X?via%3Dihub
    # recomb only: not found in native, but may be found in ecoli
    # ecoli only: not found in recomb or native mammalian proteins, only in ecoli
    known_substitutions  =                ['A to T',
                                           'A to V', 
                                           'D to G', # recomb only, # PTM
                                           'D to E',                # PTM
                                           'D to N',                # PTM
                                           'E to K',
                                           'F to I',
                                           'F to L',
                                           'G to D', # recomb only
                                           'G to E',             # ecoli only
                                           'G to R',             # ecoli only
                                           'G to S',             # ecoli only
                                           'H to Q',             # ecoli only
                                           'H to Y',
                                           'K to R', # recomb only
                                           'L to F',
                                           'L to K',              # ecoli only
                                           'L to V',              # ecoli only
                                           'M to I',
                                           'M to L',
                                           'M to T',
                                           'N to K', # recomb only
                                           'N to S', # recomb only
                                           'P to I',            # ecoli only
                                           'P to L',            # ecoli only
                                           'P to S',              # ecoli only
                                           'Q to H',              # ecoli only
                                           'R to G', # recomb only
                                           'R to K',
                                           'R to Q',              # ecoli only 
                                           'S to I', 
                                           'S to L', 
                                           'S to N', # recomb only
                                           'S to R', # recomb only
                                           'S to T',              # ecoli only
                                           'T to S', # recomb only
                                           'T to I',
                                           'T to L',
                                           'V to A', 
                                           'V to I', 
                                           'V to L', 
                                           'V to M',              # ecoli only
                                           'Y to F', # recomb only
                                           'Y to H', # recomb only
                                           'Y to N']              # ecoli] 
    
    recombinant_only_substitutions = [     'D to G', # recomb only
                                           'G to D', # recomb only
                                           'K to R', # recomb only
                                           'N to K', # recomb only
                                           'N to S', # recomb only
                                           'R to G', # recomb only
                                           'S to N', # recomb only
                                           'S to R', # recomb only
                                           'T to S', # recomb only
                                           'Y to F', # recomb only
                                           'Y to H'] # recomb only
    
    ecoli_only_substitutions = [           'G to E',             # ecoli only
                                           'G to R',             # ecoli only
                                           'G to S',             # ecoli only
                                           'H to Q',             # ecoli only
                                           'L to K',              # ecoli only
                                           'L to V',              # ecoli only
                                           'P to L',            # ecoli only
                                           'P to I',            # ecoli only
                                           'P to S',              # ecoli only
                                           'Q to H',              # ecoli only
                                           'R to Q',              # ecoli only 
                                           'S to T',              # ecoli only
                                           'V to M',              # ecoli only
                                           'Y to N']              # ecoli only
    # potential PTMs, ignoring positions
    potential_PTMs = [  
        'A to E',
        'A to G',
        'A to Q',
        'A to S',
        'D to A',
        'D to G',
        'D to E',
        'D to H',
        'D to N',
        'E to A',
        'E to D',
        'E to Q',
        'F to E',
        'F to V',
        'F to Y',
        'G to A',
        'G to D',
        'G to N',
        'H to N',
        'I to V',
        'I to N',
        'I to Q',
        'I to R',
        'L to V',
        'L to N',
        'L to Q',
        'L to R',
        'K to R',
        'K to W',
        'M to S',
        'N to G',
        'N to P',
        'P to E',
        'P to Q',
        'P to S',
        'P to T',
        'Q to A',
        'Q to N',
        'R to K',
        'S to G',
        'T to G',
        'T to A',
        'T to D',
        'T to S',
        'V to F',
        'V to M',
        'V to N',
        'V to D',
        'Y to D']
    # 
    uniq_subs['known_substitution'] = ''
    uniq_subs['potential_PTM'] = ''
    n_known = 0
    n_ptm = 0
    for i in uniq_subs.index:
        subs = uniq_subs.at[i,'origin']+' to '+uniq_subs.at[i,'destination']
        if subs in known_substitutions:
            uniq_subs.at[i,'known_substitution'] = 'native'
            n_known += 1
            if subs in recombinant_only_substitutions:
                uniq_subs.at[i,'known_substitution'] = 'recombinant_only'
            elif subs in ecoli_only_substitutions:
                uniq_subs.at[i,'known_substitution'] = 'ecoli_only'
        if subs in potential_PTMs: 
            uniq_subs.at[i,'potential_PTM'] = 'potential_PTM'
            n_ptm += 1
    print("- marked as known substitution types: "+str(n_known))
    print("- marked as potential PTM without position filter: "+str(n_ptm))
    
    # SNPs
    uniq_subs['SNV'] = 0
    n=0
    if os.path.isfile(path_to_variant_peptide):
        variant_pep = open(path_to_variant_peptide).read().replace("\n","")
        for i in uniq_subs.index:
            if uniq_subs.at[i,'modified_sequence'] in variant_pep:
                uniq_subs.at[i,'SNV'] = 1
                n=n+1
        print("- peptides marked as potential SNP peptides: "+str(n))
    else:
        print("- skip variant filter: variant peptide file not found or not set: "+path_to_variant_peptide)
    uniq_subs.to_csv(path_to_uniq_subs[:-4]+'-filtered.csv')
        
def plotting(output_dir):
    
    uniq_subs = pd.read_csv(output_dir+'/uniq_subs-filtered.csv')

    # histogram of error rate
    nbin = 30
    
    fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, tight_layout=True)
    d = uniq_subs[uniq_subs.longest_protein_status == 'normal']
    _, bins, _ = axs[0].hist(d.log10_mod_to_ref_intensity_ratio,nbin,alpha=.3,label='Normal,'+str(len(d)))
    axs[0].hist(d.log10_mod_to_ref_intensity_ratio[d.mispairing==1],bins=bins,alpha=.3,label='near-cognate,'+str(sum(d.mispairing==1)))
    axs[0].hist(d.log10_mod_to_ref_intensity_ratio[d.SNV==1],bins=bins,alpha=.3,label='SNV,'+str(sum(d.SNV==1)))
    axs[0].legend()
    
    # antibody
    d = uniq_subs[uniq_subs.longest_protein_status != 'normal']
    axs[1].hist(d.log10_mod_to_ref_intensity_ratio,bins=bins,alpha=.3,label='Antibody,'+str(len(d)))
    axs[1].hist(d.log10_mod_to_ref_intensity_ratio[d.mispairing==1],bins=bins,alpha=.3,label='near-cognate,'+str(sum(d.mispairing==1)))
    axs[1].hist(d.log10_mod_to_ref_intensity_ratio[d.SNV==1],bins=bins,alpha=.3,label='SNV,'+str(sum(d.SNV==1)))
    axs[1].legend()
    plt.savefig(output_dir+'/uniq-subs-intensity-ratio-hist.pdf')  
    
    # other plots only look at non SNPs, with log10_mod_to_ref_intensity_ratio < -1
    # in aggregate plots, can use SNPs as background to see bias in positions like 3rd codon position GU wobble pair?
    d = uniq_subs[uniq_subs.log10_mod_to_ref_intensity_ratio < -1]
    
    # plot aa to aa heatmap
    aa2aa = aa_to_aa_count_matrix(d)
    plot_matrix(np.log2(aa2aa+1),output_dir+'/uniq-subs-aa2aa-log-heatmap.pdf',10)
    aa2aa.to_csv(output_dir+'/uniq-subs-aa2aa-matrix.csv', header=True, index=True, float_format='%d')

    # plot codon to aa heatmap
    d = uniq_subs[uniq_subs.log10_mod_to_ref_intensity_ratio < -1]
    codon2aa = codon_to_aa_count_matrix(d)
    #codon2aa = codon_to_aa_count_matrix(uniq_subs,'','') # do not mask stop codons or PTMs
    #codon2aa = codon_to_aa_count_matrix(uniq_subs,'mask stop codons','mask PTM') # do not mask stop codons or PTMs
    #codon2aa = codon_to_aa_count_matrix(uniq_subs,'delete stop codons','mask PTM') # do not mask stop codons or PTMs

    plot_matrix(np.log2(codon2aa+1),output_dir+'/uniq-subs-codon2aa-log-heatmap.pdf')
    codon2aa.to_csv(output_dir+'/uniq-subs-codon2aa-matrix.csv', header=True, index=True, float_format='%d')
    
    # plot codon to aa heatmap, mispairing only
    d = uniq_subs[uniq_subs.log10_mod_to_ref_intensity_ratio < -1]
    codon2aa = codon_to_aa_count_matrix(d[d['mispairing']==1])
    plot_matrix(np.log2(codon2aa+1),output_dir+'/uniq-subs-codon2aa-mispairing-only-log-heatmap.pdf')
    codon2aa.to_csv(output_dir+'/uniq-subs-codon2aa-mispairing-only-matrix.csv', header=True, index=True, float_format='%d')

def plot_matrix(m,filename,fontsize=5):
    fig, ax = plt.subplots()
    im = ax.imshow(m,cmap='magma')

    ax.set_facecolor("gray")

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(m.columns)))
    ax.set_yticks(np.arange(len(m.index)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(m.columns,fontsize=fontsize)
    ax.set_yticklabels(m.index,fontsize=fontsize)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    '''
    for i in range(len(index)):
        for j in range(len(column)):
            text = ax.text(j, i, m[i, j],ha="center", va="center", color="w")
    '''

    #ax.set_title("Title")
    fig.tight_layout()
    #plt.show()
    plt.savefig(filename)  

    '''
    # example
    index = ["AAA", "CCC", "UUU", "GGG", "CGA", "GUG", "UGG"]
    columns = list('ABCDEFG')
    d = np.array([[0.8, 2.4, 2.5, 3.9, 0.0, 4.0, 0.0],
                            [2.4, 0.0, 4.0, 1.0, 2.7, 0.0, 0.0],
                            [1.1, 2.4, 0.8, 4.3, 1.9, 4.4, 0.0],
                            [0.6, 0.0, 0.3, 0.0, 3.1, 0.0, 0.0],
                            [0.7, 1.7, 0.6, 2.6, 2.2, 6.2, 0.0],
                            [1.3, 1.2, 0.0, 0.0, 0.0, 3.2, 5.1],
                            [0.1, 2.0, 0.0, 1.4, 0.0, 1.9, 6.3]])
    m = pd.DataFrame(data = d, index = index, columns=columns,dtype=float)
    plot_matrix(m,'test.pdf')
    '''

def codon_to_aa_count_matrix(uniq_subs):
    # mask_stop_codons: 'mask stop codons' or 'delete stop codons' or else
    # mask_potential_PTMs: 'mask PTM' or else
    
    uniq_subs = uniq_subs[pd.notnull(uniq_subs['codon'])]
    uniq_subs['codon'] = uniq_subs['codon'].map(lambda x: x.replace('T','U'))
    counts = uniq_subs[['codon','destination','origin']].groupby(['codon','destination']).count()	
    #print(counts.index)

    bases = 'UCAG'
    codons = [a+b+c for a in bases for b in bases for c in bases]
    matrix = pd.DataFrame(data = 0, index = codons, columns=list('ACDEFGHKLMNPQRSTVWY'),dtype=float)
    
    
    for codon,destination in counts.index:
        matrix.loc[codon,destination] = counts.origin[codon][destination]
    
    #print(matrix.index)
    
    # the following code will gray out stop codons or delete them
    for label in matrix.index:
        ## the following code will gray out stop codons
        if label == 'UGA' or label == 'UAG' or label == 'UAA':
            #if mask_stop_codons == 'mask stop codons':
            matrix.loc[label] = 0 #float('NaN')
            #elif mask_stop_codons == 'delete stop codons':
            #    matrix.drop(label)
        
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG' 
    codon_table = get_codon_table(codons,amino_acids)
    inverted_codon_table = get_inverted_codon_table(codons,amino_acids)
    inverted_codon_table['L'] = inverted_codon_table['L'] + inverted_codon_table['I']

    exact_PTM_spec_list =   ['A to I',
                                 'A to L',
                                 'S to T',
                                 'S to D',
                                 'S to E',
                                 'P to E',
                                 'T to E',
                                 'N to D',
                                 'N to Q',
                                 'N to C',
                                 'D to E',
                                 'Q to E',
                                 'K to C',
                                 'F to C',
                                 'F to Y']


    ## the following code will gray out potential PTM substitutions
    #print(matrix.index)
    for label in matrix.index:
        for col in matrix.columns:
            if (label in inverted_codon_table[col]) or (codon_table[label] +' to '+col in exact_PTM_spec_list):
                matrix.loc[label, col] = 0#float('NaN')
    
    matrix.rename(columns={"L": "I\nL"},inplace=True)
    
    # add aa to codons as row name
    for label in matrix.index:
        matrix.at[label,'new-index'] = label+'-'+codon_table[label]
    matrix.set_index('new-index',inplace=True)
    
    return matrix

def aa_to_aa_count_matrix(uniq_subs):
    exact_PTM_spec_list =   ['A to I',
                             'A to L',
                             'S to T',
                             'S to D',
                             'S to E',
                             'P to E',
                             'T to E',
                             'N to D',
                             'N to Q',
                             'N to C',
                             'D to E',
                             'Q to E',
                             'K to C',
                             'F to C',
                             'F to Y']
    matrix = pd.DataFrame(data = 0, index =list('ACDEFGHKLMNPQRSTVWY'), columns=list('ACDEFGHKLMNPQRSTVWY'),dtype=float)
    for i in uniq_subs.index:
        matrix.loc[uniq_subs.at[i,'origin'].replace('I','L'),uniq_subs.at[i,'destination']] += 1
    matrix.rename(columns={"L": "I\nL"},inplace=True)
    matrix.rename(index={"L": "I/L"},inplace=True)
    
    for label in matrix.index:
        for col in matrix.columns:
            if label +' to '+col in exact_PTM_spec_list:
                matrix.loc[label, col] = 0#float('NaN')
                
    return matrix

def merge(folder):
    # merge the uniq_subs file from multiple samples
    data = []
    for f in os.scandir(folder):
        if f.is_dir():
            sample = f.name
            if os.path.isfile(folder+'/'+sample+'/substitution/uniq_subs-filtered.csv'):
                uniq_subs = pd.read_csv(folder+'/'+sample+'/substitution/uniq_subs-filtered.csv')
                uniq_subs['sample'] = sample
                data.append(uniq_subs)
                print(sample)
    alldata = pd.concat(data)
    alldata.to_csv(folder+'/uniq_subs-filtered.csv')
    plotting(folder)

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

    parser.add_argument('input_dir',help='Required. Path to a folder with raw data (*.raw) (exception: see --merge)')

    parser.add_argument("--analysis", dest='analysis', action='store',default='canonical',
                         help='combination of canonical, frameshift, utr, lncrna, intron, and substitution separated by comma (no space). Default: canonical')

    parser.add_argument("--output-dir",action='store',default='NA',
                         help='Output folder name. Default: same as input')

    parser.add_argument("--transcriptome", action='store',default=transcriptome,
                         help='Path to transcriptome fasta file. See README for details')
              
    parser.add_argument("--proteome", action='store',default=proteome,
                         help='Path to proteome fasta file')
    
    parser.add_argument("--contaminant", action='store',default=contaminant,
                         help='Path to contaminant fasta file')
    
    parser.add_argument("--variant", action='store',default=variant,
                         help='Path to variant peptides created by SNV. To skip, set it NA')

    parser.add_argument("--substitution-xml", dest='template_xml_substitution', action='store',default=template_xml_substitution,
                         help='A template xml file for substitution detection')
    
    parser.add_argument("--standard-xml", dest='template_xml_standard', action='store',default=template_xml_standard,
                         help='A template xml file for standard MaxQuant search')
    
    parser.add_argument("--thread", dest='thread', action='store',default=4,type=int,
                         help='Number fo threads for MaxQuant. Default: 4')
    
    parser.add_argument("--merge", action='store_true',
                         help='merge output from multiple substitution analysis. Each folder under input dir is one sample')
              
    args = parser.parse_args()
    
    if args.merge:
        merge(args.input_dir)
        exit(0)

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

    # start the analysis pipeline
    analyses = args.analysis.lower().split(',')

    # for analysis other than 'substitution'
    for analysis in analyses:
        
        if analysis == 'substitution':
            continue
            
        # create sub-folder for specific type of analysis (canonical,frameshift, utr, lncrna, intron etc)
        output_dir = args.output_dir +'/'+analysis         
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        
        # generate custom proteome if not aligning to the canonical proteome
        path_to_custom_proteome = args.proteome
        if analysis != 'canonical':
            path_to_custom_proteome = generate_custom_proteome(analysis)

        # print to screen the parameters
        print("=========================================")
        print("analysis         : "+analysis)
        print("=========================================")
        print("- proteome       : "+path_to_custom_proteome)
        print("- variant        : "+args.variant)
        print("- template xml   : "+args.template_xml_standard)
        print("- thread         : "+str(args.thread))
        print("- output         : "+output_dir)
        print("- input          : "+args.input_dir)
        print("                   " + str(nSample) + " raw file(s)")

        # run maxquant
        maxquant_standard_search(path_to_custom_proteome,args.input_dir,output_dir,args.template_xml_standard,args.thread)
        
        if analysis != 'canonical':
            npep = filter_maxquant_result(output_dir,args.proteome,args.variant)
        
    # substitution 
    
    if 'substitution' in args.analysis.lower():
        
        args.output_dir = args.output_dir +'/substitution'

        print("Path to files:")
        print("- proteome          : "+args.proteome)
        print("- variant           : "+args.variant)
        print("- transcriptome     : "+args.transcriptome)
        print("- contaminant       : "+args.contaminant)
        print("- template xml      : "+args.template_xml_substitution)
        print("- thread            : "+str(args.thread))
        print("- output            : "+args.output_dir)
        print("- input             : "+args.input_dir)
        print("                      " + str(nSample) + " raw file(s)")

        # dependent peptide search by MaxQuant
        if os.path.isfile(args.output_dir+'/combined/txt/allPeptides.txt'):
            print("MaxQuant search results detected in the output folder "+args.output_dir)
            print("- to run MaxQuant again, please change output directory")
        else:
            print("Generate MaxQuant parameter file (xml)")
            generate_xml(args.template_xml_substitution,args.input_dir,args.output_dir,args.proteome,args.thread)
            print("MaxQuant search of dependent peptides")
            cmd = MaxQuantCmd+ " " + args.output_dir + "/mqpar.xml"
            print('- '+cmd)
            os.system(cmd)
            print("Deleting intermediate files from MaxQuant run")
            clean_up_after_maxquant_run(args.input_dir,args.output_dir)
            
        # identify potential substitutions
        if os.path.isfile(args.output_dir+'/subs.csv'):
            print("Substitution detection results found ")
            print("- to run again, please delete the file "+args.output_dir+'/subs.csv')
        else:
            print("Substitution detection and filtering")
            substitution_detection_and_filtering(args.input_dir,args.output_dir,args.transcriptome)
            #import detect
            #import quantify
            #import plot
            
        # standard maxquant search to get intensity of substitutions
        print("Second pass")
        substitution_second_pass(args.output_dir+'/subs.csv',args.input_dir,args.output_dir,args.template_xml_standard,args.thread)
        print("Substitution quantification")
        substitution_quantification_with_2nd_pass(args.output_dir+'/subs.csv',args.output_dir+'/second_pass/combined/txt/evidence.txt',args.output_dir)
        print("Additional filters")
        substitution_filter(args.output_dir+'/uniq_subs.csv',args.variant,args.proteome,args.contaminant)
        print("Plotting results")
        plotting(args.output_dir)

        