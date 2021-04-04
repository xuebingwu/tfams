# generate a fasta file for variant peptides encoded by SNPs
# the output will be used to filter peptides

# example usage using dbSNP as input

# python variant.py /media/rna/xuebing/MaxQuant/All_20180418.vcf.gz /media/rna/xuebing/MaxQuant/gencode.v37.pc_translations.fa --tag NSM

import os
import sys
import pandas as pd
from Bio import SeqIO
import argparse 
import re

path_to_vep='~/software/ensembl-vep/vep'

def find_missense_variant_position_in_protein_with_ensembl_vep(vcf_file,missense_filter=""):
    '''
    input:
        vcf_file
            - vcf file containing variants, such as those from dbSNP or called from RNA-seq/exomes. Can be *.vcf.gz or *.vcf
        missense_filter
            - keep lines containing specified tag, such as 'NSM' (missense variant) for VCF files from dbSNP database
    output:
        <vcf_file>.missense.vep.txt
            - VEP annotated output
    '''
    
    if os.path.isfile(vcf_file):
        
        # if output file already exists, do nothing
        if os.path.isfile(vcf_file+'.missense.vep.txt'):
            print("VEP output file exists: "+vcf_file+'.missense.vep.txt')
            return 0
        if os.path.isfile(vcf_file+'.missense.vep.txt.gz'):
            print("VEP output file exists: "+vcf_file+'.missense.vep.txt.gz')
            return 0       
        
        # otherwise, run VEP
        vep_input_file = vcf_file
        if missense_filter != "":
            print("Filtering vcf file")
            cmd = "cat " + vcf_file + " | grep -e '#' -e " + missense_filter + " > " + vcf_file + "." + missense_filter + ".vcf"
            if vcf_file[-3:] == '.gz':
                cmd = 'z'+cmd
            print("    " + cmd)
            os.system(cmd)
            vep_input_file = vcf_file + "." + missense_filter + ".vcf"
            if not os.path.isfile(vep_input_file):
                print(" - vcf file filtering failed: output file not found: "+vep_input_file)
        elif vcf_file[-3:] == '.gz':
            print("Unzipping vcf file")
            vep_input_file = vcf_file+'.vcf.tmp'
            cmd = 'zcat '+vcf_file+' > '+vep_input_file
            print("    " + cmd)
            os.system(cmd)
            
        # annotate variants with VEP. takes hours/days
        cmd = path_to_vep + " -i " + vep_input_file + " -o " + vcf_file + ".vep.txt --cache --protein --symbol --sift b"
        print("Running VEP, can take hours/days")
        print("    " + cmd)
        os.system(cmd)
        
        # filter to only keep vep annotatted missense variants
        cmd = 'cat '+ vcf_file + '.vep.txt | grep missense > ' + vcf_file + '.missense.vep.txt'
        print("Saving missense variants")
        print("    " + cmd)
        os.system(cmd)
        
        # compress some intermediate files that might be useful   
        if os.path.isfile(vcf_file + ".vep.txt"):
            os.system('gzip '+ vcf_file + ".vep.txt" )
        if os.path.isfile(vcf_file + "." + missense_filter + ".vcf"):
            os.system('gzip '+vcf_file + "." + missense_filter + ".vcf")
        
        # deleting some intermediate files
        if os.path.isfile(vcf_file + '.vep.txt_summary.html'):
            os.system('rm '+vcf_file + '.vep.txt_summary.html')
        if os.path.isfile(vcf_file+'.vcf.tmp'):
            os.system('rm '+vcf_file+'.vcf.tmp')
    else:
        print("ERROR: vcf file not exist: "+vcf_file)
        exit(1)

def extract_amino_acid_sequence_flanking_missense_variants(vep_missense_file,proteome_file,output_file,flanking_len=50):
    '''
    input:
        vep_output_file
            - output from the snp2pep.sh script, i.e. missense variants from VEP annotated VCF
        proteome_file
            - amino acid sequences of all proteins
            - can be downloaded from GENCODE (e.g. gencode.v37.pc_translations.fa)
            - header needs to be separted by '|' and has ENSEMBL protein id in field 1, i.e. >ENSP00000478421.2|(other-info)
            - the one from uniprot will not work.
        output_file
            - output file path and name
        flanking_len
            - max lenght of amino acids on each side of the variant to be included in the output
            - should be longer than the min peptide length in MS, ideally the max peptide length
            - recommended: 50
    output:
        a fasta file
    '''
    
    input_file = vep_missense_file
    if not os.path.isfile(input_file):
        if os.path.isfile(input_file+'.gz'):
            os.system('zcat '+input_file+'.gz > ' + input_file + '.tmp ')
            input_file = input_file + '.tmp'
        else:
            print("ERROR: VEP output file not found: "+input_file)
            exit(1)
    
    print("Extracting flanking sequences")
    
    # load protein sequences into a dictionary
    Seqs = {}
    for record in SeqIO.parse(open(proteome_file,'r'),'fasta'):
        descriptions = record.description.split('|') 
        # >ENSP00000478421.2|ENST00000616016.5|ENSG00000187634.13|OTTHUMG00000040719.11|OTTHUMT00000316521.3|SAMD11-209|SAMD11|844
        Seqs[descriptions[0].split('.')[0]] = record.seq.upper()   
    print(str(len(Seqs))+' protein sequences loaded')

    # get protein sequences around the variants
    # rs140739101     1:69428 G       ENSG00000186092 ENST00000335137 Transcript      missense_variant        374     338     113     F/C     tTt/tGt -       IMPACT=MODERATE;STRAND=1;ENSP=ENSP00000334393;SIFT=deleterious(0.04)

    out = open(output_file,'w')
    n_variant = 0      # total number of variants processed
    n_variant_protein_not_found = 0 # total number of variants skipped due to not-found protein
    n_variant_wrong_reference  = 0 # referene allele wrong
    n_variant_no_sift = 0
    proteins_not_found = {}
    with open(input_file) as f:
        for line in f:
            n_variant += 1
            # get variant info: protein, position, and mutations
            flds = line.strip().split('\t')
            variant = flds[0]
            protein = re.findall(r'ENSP=(.*?);',line)
            gene = re.findall(r'SYMBOL=(.*?);',line)
            # note that for NMD substrates there won't be SIFT, so no ; after ENSP ID. Ignore those
            if (len(protein) == 0) or (len(gene) == 0):
                n_variant_no_sift += 1
                continue
            else:
                protein = protein[0]
                gene = gene[0]
            
            position = int(flds[9])-1
            source,dest = flds[10].split('/') # will there be more than 2 alleles
            SIFT = line.strip().split(';')[-1]
            
            # protein not found
            if not (protein in Seqs):
                proteins_not_found[protein] = 1
                n_variant_protein_not_found += 1
                continue

            # check if the sequence match, and then make the mutant sequence 
            if Seqs[protein][position] == source:
                start = max(0,position-flanking_len)
                end = position+flanking_len
                seq = Seqs[protein][start:position]+dest+Seqs[protein][position+1:end]
                out.write('>'+gene+'|'+protein+'|'+variant+'|'+str(start)+'|'+source+'|'+dest+'|'+SIFT+'\n'+str(seq)+'\n')
            else:
                print("incorrect position/sequence for variant "+variant)
                n_variant_wrong_reference += 1
    print(str(n_variant) + " missense variants processed")
    print(str(len(proteins_not_found)) + " proteins not found")
    print(str(n_variant_protein_not_found) + " variants in proteins not found")
    print(str(n_variant_no_sift) + " variants without SIFT predictions discarded")
    if n_variant_wrong_reference>0:
        print(str(n_variant_wrong_reference) + " variants with wrong reference allele discarded ")
    print(str(n_variant - n_variant_protein_not_found - n_variant_wrong_reference - n_variant_no_sift) + " variants saved in "+output_file)
    
    # delete unzipped input file
    if input_file[-4:] == '.tmp':
        os.system('rm '+input_file)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    parser.add_argument('variant',help='Required. Path to a variant file in VCF format')
    parser.add_argument('proteome',help='Required. Path to GENCODE protein sequence fasta file')
    parser.add_argument("--output",action='store',default='NA',
                         help='Output file name. Default: [input]-snp2pep.fa')
    parser.add_argument("--length", action='store',default=50,
                         help='Max number of amino acids flanking the variant to output. Max total length is 2N+1')
    parser.add_argument("--tag",action='store',default='',
                         help='A tag for filtering input variants, such as NSM (missense variant) for VCF files from dbSNP')
    
    args = parser.parse_args()
            
    if args.output == 'NA':
        args.output = args.variant+'-snp2pep.fa'
    
    find_missense_variant_position_in_protein_with_ensembl_vep(args.variant,missense_filter=args.tag)
    
    extract_amino_acid_sequence_flanking_missense_variants(args.variant+'.missense.vep.txt',args.proteome,args.output,args.length)
    
    # compress vep output
    if os.path.isfile(args.variant+'.missense.vep.txt'):
        os.system('gzip '+ args.variant+'.missense.vep.txt' )