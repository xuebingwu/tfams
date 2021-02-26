# TFAMS
## Translation Fidelity Analysis with Mass Spectrometry Data

This is a pipeline for quantifying translation fidelity with mass spectrometry data. The pipeline provides a wrapper for streamlining MaxQuant dependent peptide search in Linux and subsequent detection and quantification of amino acid substitutions, focusing on near-cognate substitutions likely caused by misdecoding. The detection/quantification part is modified from scripts written by Ernest Mordret (https://github.com/ernestmordret/substitutions/).

## Requirement / setup

### Python packages

The pipeline is written in Python 3 and requires the following packages:

```numpy, scipy, pandas, matplotlib (for plotting), biopython```

### MaxQuant

The pipeline uses MaxQuant to search for potential amino acid substitutions. To allow MaxQuant to run in Linux commandline, follow the instructions here to download ```MaxQuant``` and install ```.NET CORE```: http://coxdocs.org/doku.php?id=maxquant:common:download_and_installation

Once installed, please update the path to ```MaxQuantCmd.exe``` in the beginning of the script ```tfams.py```:

```
MaxQuantCmd="dotnet MaxQuant/bin/MaxQuantCmd.exe"
``` 

### Template parameter file (xml) for MaxQuant
Two template xml files are provided in the folder ```./reference```. Typically no need to change it.

### Reference proteome and transcriptome

Proteome: The amino acid sequences for all proteins in fasta format. Used by MaxQuant and can be downloaded from UniProt: https://www.uniprot.org/proteomes/ or GENCODE: https://www.gencodegenes.org/

Transcriptome: The coding sequence (CDS) of all mRNAs in fasta format. Can be downloaded from ENSEMBL: http://ftp.ensembl.org/pub/release-103/fasta/

To setup the default paths to the proteome and transcriptome, edit the following line at the beginning of the script ```tfams.py```:

```
proteome='./reference/human.protein.fa'
transcriptome='./reference/human.CDS.fa'
template_xml='./template_xml/mqpar-human.xml'
```
## Usage

```
usage: tfams.py [-h] --input INPUT_DIR [--output OUTPUT_DIR] [--proteome PROTEOME] [--transcriptome TRANSCRIPTOME] [--xml TEMPLATE_XML]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT_DIR     Required. Path to a folder with raw data (*.raw)
  --output OUTPUT_DIR   Output folder name. Default: same as input
  --proteome PROTEOME   Path to proteome fasta file
  --transcriptome TRANSCRIPTOME
                        Path to transcriptome (CDS only) fasta file
  --xml TEMPLATE_XML    A template xml file
``` 

### Example

Analyze human data with default settings: 

```
python tfams.py --input raw_file_folder
``` 

