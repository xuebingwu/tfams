# TFAMS
## Translation Fidelity Analysis with Mass Spectrometry Data

TFAMS is a pipeline for detecting translation errors from mass spectrometry data. Types of errors include: 1) amino acid substitutions likely caused by tRNA mispairing, 2) frameshifting, 3) translation in noncoding sequences, including 5' UTRs, 3' UTRs, introns, and long noncoding RNAs (lncRNAs). The amino acid substitution part is modified from scripts written by Ernest Mordret (https://github.com/ernestmordret/substitutions/).

## Requirement / setup

The pipeline is written in Python 3 and has only been tested in Linux (Ubuntu 20.04.2 LTS).

### Python packages

`numpy, scipy, pandas, biopython, matplotlib (for plotting)`

### MaxQuant

The pipeline uses MaxQuant to search for peptide search. To allow MaxQuant to run in Linux commandline, follow the instructions here to download `MaxQuant` and install `.NET CORE`: http://coxdocs.org/doku.php?id=maxquant:common:download_and_installation

Once installed, please update the path to `MaxQuantCmd.exe` in the beginning of the script `substitution.py` and `noncanonical_translation.py`:

```python
MaxQuantCmd="dotnet MaxQuant/bin/MaxQuantCmd.exe"
``` 

### Template parameter file (xml) for MaxQuant
Template xml files are provided in the folder `./template_xml`. Typically no need to change.

### Reference proteome and transcriptome

Proteome: The amino acid sequences for all proteins in fasta format. Used by MaxQuant and can be downloaded from UniProt: https://www.uniprot.org/proteomes/ or GENCODE: https://www.gencodegenes.org/

Transcriptome: For substution or frameshift detection the transcriptome is the coding sequence (CDS only, no UTRs) of all mRNAs in fasta format. Can be downloaded from ENSEMBL: http://ftp.ensembl.org/pub/release-103/fasta/. For others see details below.

To setup the default paths to the proteome and transcriptome files, edit the following lines at the beginning of the script `substitution.py`:

```python
transcriptome='./reference/human.CDS.fa'
proteome='./reference/human.protein.fa'
```

and `noncanonical_translation.py`:

```python
transcriptome_frameshift='./reference/human.CDS.fa'    # for frameshift, same as substitution
transcriptome_lncrna='./reference/human.lncRNA.fa'     # for lncRNA analysis, downloaded from GENCODE, lncRNA sequence
transcriptome_mrna='./reference/human.mRNA.fa'         # for UTR analysis, downloaded from GENCODE, protein-coding transcript sequence
transcriptome_intron='./reference/human.intron.fa'     # for intron analysis, downloaded from UCSC Table browser, gencode.v32, +9nt flanking sequence
proteome='./reference/human.protein.fa'
```

## Usage

### substitution detection

Analyze human data with default settings: 

```sh
python substitution.py raw_file_folder
``` 

Detailed usage:

```sh
usage: substitution.py [-h] [--output-dir OUTPUT_DIR] [--proteome PROTEOME] [--transcriptome TRANSCRIPTOME] [--template-xml TEMPLATE_XML]
                       input_dir

positional arguments:
  input_dir             Required. Path to a folder with raw data (*.raw)

optional arguments:
  -h, --help            show this help message and exit
  --output-dir OUTPUT_DIR
                        Output folder name. Default: same as input
  --proteome PROTEOME   Path to proteome fasta file
  --transcriptome TRANSCRIPTOME
                        Path to transcriptome (CDS only) fasta file
  --template-xml TEMPLATE_XML
                        A template xml file for substitution detection (provided in ./reference)
``` 
### Noncanonical translation 

Detect canonical peptides:

```sh
python noncanonical_translation.py raw_file_folder
```

Detect both canonical and all noncanical peptides:

```sh
python noncanonical_translation.py raw_file_folder --analysis frameshift,utr,lncrna,intron
```

Detailed usage:

```sh
python noncanonical_translation.py -h
usage: noncanonical_translation.py [-h] [--analysis ANALYSIS] [--output-dir OUTPUT_DIR] [--transcriptome TRANSCRIPTOME]
                                   [--proteome PROTEOME] [--template-xml TEMPLATE_XML]
                                   input_dir

positional arguments:
  input_dir             Required. Path to a folder with raw data (*.raw)

optional arguments:
  -h, --help            show this help message and exit
  --analysis ANALYSIS   combination of frameshift, utr, lncrna, or intron separated by comma (no space)
  --output-dir OUTPUT_DIR
                        Output folder name. Default: same as input
  --transcriptome TRANSCRIPTOME
                        Path to transcriptome fasta file. See README for details
  --proteome PROTEOME   Path to proteome fasta file
  --template-xml TEMPLATE_XML
                        A template xml file
``` 
