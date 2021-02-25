# TFAMS
## Translation Fidelity Analysis with Mass Spectrometry Data

This is a pipeline for quantifying translation fidelity with mass spectrometry data. The pipeline provides a wrapper for streamlining MaxQuant dependent peptide search in Linux and subsequent detection and quantification of amino acid substitutions. The latter part is modified from scripts written by Ernest Mordret (https://github.com/ernestmordret/substitutions/), which was written in Python 2 and used for analyzing E.coli and yeast data.

## Requirement

### Python packages

The pipeline is written in Python 3 and requires the following packages:

```numpy, scipy, pandas, matplotlib (for plotting), and biopython```

### MaxQuant

To run MaxQuant from commandline, follow the instructions here to download ```MaxQuant``` and install ```.NET CORE```: http://coxdocs.org/doku.php?id=maxquant:common:download_and_installation

Once installed, please update the following line in the script ```tfams.py``` with path to ```MaxQuantCmd.exe```:

```
MaxQuantCmd="dotnet MaxQuant/bin/MaxQuantCmd.exe"
``` 

### Reference proteome and transcriptome

Proteome: The amino acid sequences for all proteins in fasta format. Used by MaxQuant and usually downloaded from UniProt: https://www.uniprot.org/proteomes/

Transcriptome: The coding sequence (CDS) of all mRNAs in fasta format. Can be downloaded from ENSEMBL: http://ftp.ensembl.org/pub/release-103/fasta/


## Usage

```
usage: tfams.py [-h] [--input INPUT_DIR] [--output OUTPUT_DIR] [--proteome PROTEOME] [--transcriptome TRANSCRIPTOME] [--xml TEMPLATE_XML]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT_DIR     Required. Path to a folder with raw data (*.raw)
  --output OUTPUT_DIR   Output folder name. Default: same as input
  --proteome PROTEOME   path to proteome fasta file (default: ./reference/human.protein.fa)
  --transcriptome TRANSCRIPTOME
                        path to transcriptome fasta file (default: ./reference/human.CDS.fa)
  --xml TEMPLATE_XML    A template xml file (default: ./template_xml/mqpar-human.xml)
``` 

### Example

Analyze human data with default settings: 

```
python tfams.py --input raw_file_folder
``` 

