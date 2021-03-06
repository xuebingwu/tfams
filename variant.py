# generate a fasta file for variant peptides encoded by SNPs

import os
import sys
import pandas as pd
from Bio import SeqIO
import argparse


def extract_amino_acid_sequence_flanking_missense_variants(variant_file,proteome_file,output_file,flanking_len):
    '''
    input:
        variant_file
            - output from the snp2pep.sh script, i.e. missense variants from VEP annotated VCF
        proteome_file
            - amino acid sequences of all proteins
            - can be downloaded from GENCODE (e.g. gencode.v37.pc_translations.fa)
            - header needs to be separted by '|' and has ENSEMBL protein id in field 1, i.e. >ENSP00000478421.2|(other-info)
        output_file
            - output file path and name
        flanking_len
            - max lenght of amino acids on each side of the variant to be included in the output
            - should be longer than the min peptide length in MS, ideally the max peptide length
            - recommended: 50
    output:
        a fasta file
    '''
    # load protein sequences into a dictionary
    Seqs = {}
    for record in SeqIO.parse(open(proteome_file,'r'),'fasta'):
        descriptions = record.description.split('|') 
        # >ENSP00000478421.2|ENST00000616016.5|ENSG00000187634.13|OTTHUMG00000040719.11|OTTHUMT00000316521.3|SAMD11-209|SAMD11|844
        Seqs[descriptions[0].split('.')[0]] = record.seq.upper()   
    print(str(len(Seqs))+' protein sequences loaded')

    # get protein sequences around the variants
    # rs140739101     1:69428 G       ENSG00000186092 ENST00000335137 Transcript      missense_variant        374     338     113     F/C     tTt/tGt -       IMPACT=MODERATE;STRAND=1;ENSP=ENSP00000334393

    out = open(output_file,'w')
    n_variant = 0      # total number of variants processed
    n_variant_protein_not_found = 0 # total number of variants skipped due to not-found protein
    n_variant_wrong_reference  = 0 # referene allele wrong
    proteins_not_found = {}
    with open(variant_file) as f:
        for line in f:
            n_variant += 1
            # get variant info: protein, position, and mutations
            flds = line.strip().split('\t')
            variant = flds[0]
            protein = flds[13].split('=')[-1]
            position = int(flds[9])-1
            source,dest = flds[10].split('/') # will there be more than 2 alleles

            # protein not found
            if not protein in Seqs:
                proteins_not_found[protein] = 1
                n_variant_protein_not_found += 1
                continue

            # check if the sequence match, and then make the mutant sequence 
            if Seqs[protein][position] == source:
                start = max(0,position-flanking_len)
                end = position+flanking_len
                seq = Seqs[protein][start:position]+dest+Seqs[protein][position+1:end]
                out.write('>'+protein+'|'+variant+'|'+str(start)+'|'+source+'|'+dest+'\n'+str(seq)+'\n')
            else:
                print("incorrect position/sequence for variant "+variant)
                n_variant_wrong_reference += 1
    print(str(n_variant) + " missense variants processed")
    print(str(len(proteins_not_found)) + " proteins not found")
    print(str(n_variant_protein_not_found) + " variants in proteins not found")
    if n_variant_wrong_reference>0:
        print(str(n_variant_wrong_reference) + " variants with wrong reference allele discarded ")
    print(str(n_variant - n_variant_protein_not_found - n_variant_wrong_reference) + " variants saved in "+output_file)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    parser.add_argument('variant',help='Required. Path to input, *.vcf.gz or VEP annotated missense variants')
    parser.add_argument('proteome',help='Required. Path to a folder with raw data (*.raw)')

    parser.add_argument("--output",action='store',default='NA',
                         help='Output file name. Default: [input]-snp2pep.fa')

    parser.add_argument("--gff3", action='store',default='NA',
                         help='Path to gff3 file. Required if input is *.vcf.gz')
    
    parser.add_argument("--length", action='store',default=50,
                         help='Number of amino acids flanking the variant to output. Total length is 2N+1')
              
    args = parser.parse_args()
        
    if args.variant[-3:] == '.gz':
        if args.gff3 == 'NA':
            print("--gff3 is required if the variant is *.vcf.gz")
            exit()
        else:
            print("Run snp2pep.sh")
            #os.system('./snp2pep.sh '+args.variant+' '+args.gff3)
            args.variant = args.variant+'.missense_variant'
            
    if args.output == 'NA':
        args.output = args.variant+'-snp2pep.fa'
        
    extract_amino_acid_sequence_flanking_missense_variants(args.variant,args.proteome,args.output,args.length)