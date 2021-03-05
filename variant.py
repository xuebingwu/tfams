# generate variant peptide database

import os
import pandas as pd
from Bio import SeqIO

#wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
#wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
    
path_to_variants = 'common_all_20180418.vcf'
path_to_vep = '~/software/ensembl-vep/vep'
path_to_proteome = 'gencode.v37.pc_translations.fa'
path_to_output = 'common_all_20180418.vcf.missense_variant-pep.fa'

flanking_len = 50 # number of amino acids from each side of the variant

# get SNV within genes
if not os.path.isfile(path_to_variants+'.gene.snv'):
    os.system('cat '+path_to_variants+' | grep GENEINFO | grep "VC=SNV" > '+path_to_variants+'.gene.snv')

# annotate variants, takes very long time to run
if not os.path.isfile(path_to_variants+'.gene.snv.annotated'):
    os.system(path_to_vep+' -i '+path_to_variants+'.gene.snv'+' -o '+path_to_variants+'.gene.snv.annotated --cache --protein')

# missense variants
if not os.path.isfile(path_to_variants+'.missense_variant'):
    os.system('cat '+path_to_variants+'.gene.snv.annotated' + ' | grep missense_variant > ' + path_to_variants + '.missense_variant')

# load protein sequences into a dictionary
Seqs = {}
for record in SeqIO.parse(open(path_to_proteome,'r'),'fasta'):
    descriptions = record.description.split('|') 
    # >ENSP00000478421.2|ENST00000616016.5|ENSG00000187634.13|OTTHUMG00000040719.11|OTTHUMT00000316521.3|SAMD11-209|SAMD11|844
    Seqs[descriptions[0].split('.')[0]] = record.seq.upper()   
print(str(len(Seqs))+' protein sequences loaded')

# get protein sequences around the variants
# rs140739101     1:69428 G       ENSG00000186092 ENST00000335137 Transcript      missense_variant        374     338     113     F/C     tTt/tGt -       IMPACT=MODERATE;STRAND=1;ENSP=ENSP00000334393

out = open(path_to_output,'w')
with open(path_to_variants + '.missense_variant') as f:
    for line in f:
        # get variant info: protein, position, and mutations
        flds = line.strip().split('\t')
        variant = flds[0]
        protein = flds[13].split('=')[-1]
        position = int(flds[9])-1
        source,dest = flds[10].split('/') # will there be more than 2 alleles
        
        if not protein in Seqs:
            print(protein + ' not found')
            continue
        
        # check if the sequence match, and then make the mutant sequence 
        if Seqs[protein][position] == source:
            start = max(0,position-flanking_len)
            end = position+flanking_len
            seq = Seqs[protein][start:position]+dest+Seqs[protein][position+1:end]
            out.write('>'+protein+'|'+variant+'|'+str(start)+'|'+source+'|'+dest+'\n'+str(seq)+'\n')
        else:
            print("incorrect position/sequence for variant "+variant)
            
            
