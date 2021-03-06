# A shell script to generate a list of variant peptides centered on SNPs
# Xuebing Wu

# Usage:
#        snp2pep.sh All_20180418.vcf.gz gencode.v37.annotation.gff3
# 

# before your run this script: 
# make sure bedtools and VEP (with cache) installed and update the path here:
path_to_bedtools='bedtools'
path_to_vep='~/software/ensembl-vep/vep'

if [ "$#" -ne 2 ]; then
    echo "please provide path to both vcf.gz and gff3 files"
    echo "snp2pep.sh dbSNP.vcf.gz gencode.gff3"
    exit
fi

# input 1: compressed VCF file downloaded from dbSNP VCF file:
# wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
path_to_vcf=$1 #"All_20180418.vcf.gz" # change the name if updated

# input 2: gene annotations (gff3 file) from GENCODE
path_to_gff=$2 #'gencode.v37.annotation.gff3'

## start of the analysis

# coordinates of exons 
more $path_to_gff | grep exon > exons.gff3

# first filter of SNPs: only keep SNV found in genes, excluding indels and larger substitutions
# also convert chromosome number format: 1 -> chr1 so that it's compatible with the gff
zcat $path_to_vcf | grep "#" > $path_to_vcf.gene.snv # get the header first 
zcat $path_to_vcf | grep -v "#" | grep GENEINFO | grep "VC=SNV" | awk '{print "chr"$_}'  >> $path_to_vcf.gene.snv

# second filter: exonic snv. takes about 40min
$path_to_bedtools intersect -a $path_to_vcf.gene.snv -b exons.gff3 -u > $path_to_vcf.exonic.snv

# annotate variants with VEP. takes hours
$path_to_vep -i $path_to_vcf.exonic.snv -o $path_to_vcf.annotated.exonic.snv --cache --protein

# get missense_variant
cat $path_to_vcf.annotated.exonic.snv | grep missense_variant > $path_to_vcf.missense_variant

# remove intermediate files
if [ -f "$path_to_vcf.missense_variant" ];then
    rm exons.gff3 $path_to_vcf.gene.snv $path_to_vcf.exonic.snv $path_to_vcf.annotated.exonic.snv
fi