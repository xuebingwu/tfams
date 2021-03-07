# TFAMS
## Translation Fidelity Analysis with Mass Spectrometry Data

TFAMS is a pipeline for detecting translation errors from mass spectrometry data. Types of errors include: 1) amino acid substitutions likely caused by tRNA mispairing, 2) frameshifting, 3) translation in noncoding sequences, including 5' UTRs, 3' UTRs, introns, and long noncoding RNAs (lncRNAs). The amino acid substitution part is modified from scripts written by Ernest Mordret (https://github.com/ernestmordret/substitutions/).

## Requirement / setup

The pipeline is written in Python 3 and has only been tested in Linux (Ubuntu 20.04.2 LTS).

#### Python packages

`numpy, scipy, pandas, biopython, matplotlib (for plotting)`

#### MaxQuant

The pipeline uses MaxQuant to search for peptide search. To allow MaxQuant to run in Linux commandline, follow the instructions here to download `MaxQuant` and install `.NET CORE`: http://coxdocs.org/doku.php?id=maxquant:common:download_and_installation

Once installed, please update the path to `MaxQuantCmd.exe` in the beginning of the script `tfams.py`:

```python
MaxQuantCmd="dotnet ~/software/MaxQuant/bin/MaxQuantCmd.exe"
``` 

#### Reference proteome and transcriptome

Proteome: The amino acid sequences for all proteins in fasta format. Used by MaxQuant and can be downloaded from UniProt: https://www.uniprot.org/proteomes/ or GENCODE: https://www.gencodegenes.org/

Transcriptome: For substution or frameshift detection the transcriptome is the coding sequence (CDS only, no UTRs) of all mRNAs in fasta format. Can be downloaded from ENSEMBL: http://ftp.ensembl.org/pub/release-103/fasta/. For others see details below.

To setup the default paths to those files, edit the following lines at the beginning of the script `tfams.py`:

```python
# default of --proteome
proteome='./reference/human.protein.fa'                # amino acid sequences of all proteins

# default of --transcriptome
transcriptome='./reference/human.CDS.fa'               # for substitution, only CDS of mRNAs

# default for other transcriptomes used to generate noncanonical peptide databases 
transcriptome_frameshift='./reference/human.CDS.fa'    # for frameshift analysis, same as substitution
transcriptome_lncrna='./reference/human.lncRNA.fa'     # for lncRNA analysis, downloaded from GENCODE, lncRNA sequence
transcriptome_mrna='./reference/human.mRNA.fa'         # for UTR analysis, downloaded from GENCODE, protein-coding transcript sequence
transcriptome_intron='./reference/human.intron.fa'     # for intron analysis, downloaded from UCSC Table browser, gencode.v32, +9nt flanking sequence
```

#### Template parameter file (xml) for MaxQuant (provided)
Template xml files with MaxQuant parameters are provided in the folder `./template_xml`. Typically there is no need to change these files. The path to the default xml files can be edited in the script `tfams.py`:

```python
# default of --standard-xml
template_xml_standard='./template_xml/mqpar-standard.xml'         # Standard MaxQuant search parameters, no dependent peptide search

# default of --substitution-xml
template_xml_substitution='./template_xml/mqpar-substitution.xml' # MaxQuant parameters with dependent peptide search and match between runs
```


#### SNP peptides (optional)

Heterozygous SNPs can lead to variant peptides that will be called substitutions. To mark/remove those SNP peptides, we search against a list of variant peptides encoded by SNPs. See the script `variant.py` for how to generate such a file. Once generated, please update the default path in the script `tfams.py`:

```python
# default of --variant
variant='./reference/human.variant.fa'                 # peptides encoded by SNV
```

## Example usage

Detect canonical peptides with default settings (path to transcriptome/proteome files listed above):

```sh
python tfams.py raw_file_folder
```

Detect both canonical and all noncanonical peptides (except substitutions):

```sh
python tfams.py raw_file_folder --analysis frameshift,utr,lncrna,intron
```

Detect substitutions:

```sh
python tfams.py raw_file_folder --analysis substitution
```

## Detailed usage:

```
usage: tfams.py [-h] [--analysis ANALYSIS] [--output-dir OUTPUT_DIR] [--transcriptome TRANSCRIPTOME]
                [--proteome PROTEOME] [--variant VARIANT] [--substitution-xml TEMPLATE_XML_SUBSTITUTION]
                [--standard-xml TEMPLATE_XML_STANDARD]
                input_dir

positional arguments:
  input_dir             Required. Path to a folder with raw data (*.raw)

optional arguments:
  -h, --help            show this help message and exit
  --analysis ANALYSIS   combination of canonical, frameshift, utr, lncrna, intron, and substitution
                        separated by comma (no space). Default: canonical
  --output-dir OUTPUT_DIR
                        Output folder name. Default: same as input
  --transcriptome TRANSCRIPTOME
                        Path to transcriptome fasta file. See README for details
  --proteome PROTEOME   Path to proteome fasta file
  --variant VARIANT     Path to variant peptides caused by SNV
  --substitution-xml TEMPLATE_XML_SUBSTITUTION
                        A template xml file for substitution detection
  --standard-xml TEMPLATE_XML_STANDARD
                        A template xml file for standard MaxQuant search
``` 
