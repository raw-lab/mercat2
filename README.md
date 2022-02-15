# MerCat2: python code for versatile k-mer counter and diversity estimator for database independent property analysis (DIPA) for multi-omic analysis

![GitHub Logo](doc/mercat_workflow.jpg)

## Dependencies

MerCat runs on python3 up to version 3.9. Some of its dependencies do not support 3.10 yet.

### Standard install 

```conda create -n 
mercat2 -c bioconda fastqc fastp prodigal python=3.9
```

```
pip install mercat2
```

### external dependencies

MerCat2 can run without external dependencies based on the options used.  

Required dependencies:

- When a raw read .fastq file is given
  - fastqc [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/]
  - fastp [https://github.com/OpenGene/fastp]

- when the -prod option is used
  - prodigal [https://github.com/hyattpd/Prodigal]

These are available through BioConda.

```bash
conda install -c bioconda fastqc fastp prodigal
```

## Installing MerCat2

### Anaconda Installer

- Available via BioConda: Enable BioConda repo and run `conda install mercat2`

```bash
conda install -c bioconda mercat2
```

### PIP Installer

BioConda dependencies are not automatically installed when using pip.

```bash
pip install mercat2
```

### Source Installer

- Clone mercat2 from github

```bash
git clone https://github.com/raw-lab/mercat2.git
```

- Run install_mercat2.py to install all required dependencies

## Usage

- -i I path to input file
- -f F path to folder containing input files
- -k K k-mer length
- -n N no of cores [default = all]
- -c C minimum k-mer count [default = 10]
- -prod run prodigal on nucleotide assembled contigs
  - Must be one of ['.fa', '.fna', '.ffn', '.fasta', 'fastq']
- -s S split files into chunks of S size, in MB (default is 100MB)
- -o output folder (default is 'mercat_results' in the current working directory)
- -h, --help show this help message

Mercat assumes the input file format based on the extension provided

- raw fastq file: ['.fastq']
- nucleotide fasta: ['.fa', '.fna', '.ffn', '.fasta']
- amino acid fasta: ['.faa']

## Usage examples

### Run mercat2 on a protein file (protein fasta - '.faa')

```bash
mercat2-pipeline.py -i file-name.faa -k 3 -c 10
```

### Run mercat2 on a nucleotide file (nucleotide fasta - '.fa', '.fna', '.ffn', '.fasta')

```bash
mercat2-pipeline.py -i file-name.fna -k 3 -n 8 -c 10
```

### Run mercat2 on a nucleotide file raw data (nucleotide fastq - '.fastq')

```bash
mercat2-pipeline.py -i file-name.fastq -k 3 -n 8 -c 10
```

### Run on many samples within a folder

```bash
mercat2-pipeline.py -f /path/to/input-folder -k 3 -n 8 -c 10
```

### Run on sample with prodigal option (raw reads or nucleotide contigs - '.fa', '.fna', '.ffn', '.fasta', '.fastq')

```bash
mercat2-pipeline.py -i /path/to/input-folder -k 3 -n 8 -c 10 -prod
```

- the prodigal option runs the k-mer counter on both contigs and produced amino acids

## Outputs

- Results are stored in the output folder (default 'mercat_results' of the current working directory)
  - the 'plots' folder contains an html report with interactive plotly figures
    - If at least 3 samples are provided a PCA plot will be included in the html report
  - the 'tsv' folder contains stats tables in tab separated format
    - if protein files are given, or the -prod option, a .tsv file is created for each sample containing k-mer count, pI, Molecular Weight, and Hydrophobicity metrics
    - if nucleotide files are given a .tsv file is created for each sample containing k-mer count and GC content
  - if .fastq raw reads files are used, a 'clean' folder is created with the clean fasta file.
  - if the -prod option is used, a 'prodigal' folder is created with the amino acid .faa and .gff files

![GitHub Logo](doc/PCA.png)

## Citing Mercat

If you are publishing results obtained using MerCat2, please cite:

### CONTACT

Please send all queries to Jose Luis Figueroa III [jlfiguer@uncc.edu](mailto:jlfiguer@uncc.edu)  
Dr. Richard Allen White III [rwhit101@uncc.edu](mailto:rwhit101@uncc.edu)  
Or [open an issue](https://github.com/raw-lab/mercat2/issues)
