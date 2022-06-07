# MerCat2: python code for versatile k-mer counter and diversity estimator for database independent property analysis (DIPA) for multi-omic analysis

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/mercat2/README.html)
  
![GitHub Logo](https://github.com/raw-lab/mercat2/raw/master/doc/mercat_workflow.jpg)]

## Installing MerCat2

### Option 1: Anaconda Installer

- Available via Bioconda:

```bash
conda create -n mercat2 -c conda-forge -c bioconda mercat2
conda activate mercat2
```

### Option 2: PIP Installer

- Dependencies are not automatically installed when using pip.

```bash
pip install mercat2
```

### Option 3: Source Installer

- Clone mercat2 from github
- Run install_mercat2.py to install all required dependencies

```bash
git clone https://github.com/raw-lab/mercat2.git
cd mercat2
python install_mercat2.py
```

## Dependencies

MerCat2 runs on python 3.6 up to version 3.9. Some of its dependencies do not support 3.10 yet.

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
mercat2.py -i file-name.faa -k 3 -c 10
```

### Run mercat2 on a nucleotide file (nucleotide fasta - '.fa', '.fna', '.ffn', '.fasta')

```bash
mercat2.py -i file-name.fna -k 3 -n 8 -c 10
```

### Run mercat2 on a nucleotide file raw data (nucleotide fastq - '.fastq')

```bash
mercat2.py -i file-name.fastq -k 3 -n 8 -c 10
```

### Run on many samples within a folder

```bash
mercat2.py -f /path/to/input-folder -k 3 -n 8 -c 10
```

### Run on sample with prodigal option (raw reads or nucleotide contigs - '.fa', '.fna', '.ffn', '.fasta', '.fastq')

```bash
mercat2.py -i /path/to/input-folder -k 3 -n 8 -c 10 -prod
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

![GitHub Logo](https://github.com/raw-lab/mercat2/raw/master/doc/PCA.png)

## Notes on memory usage and speed

MerCat2 uses a substantial amount of memory when the k-mer 4 is high.  
Running MerCat2 on a personal computer using a k-mer length of ~4 should be OK. Total memory usage can be reduced using the Chunker feature (-s option), but keep in mind that in testing when the chunk size is too small (1MB) some of the least significant k-mers will get lost. This does not seem to affect the overall results, but it is something to keep in mind. Using the chunker and reducing the number of CPUs available (-n option) can help reduce memory requirements.  
  
The speed of MerCat2 can be increased when more memory or computer nodes are available on a cluster and using a chunk size of about 100Mb.

## Citing Mercat

If you are publishing results obtained using MerCat2, please cite:

### CONTACT

Please send all queries to Jose Luis Figueroa III [jlfiguer@uncc.edu](mailto:jlfiguer@uncc.edu)  
Dr. Richard Allen White III [rwhit101@uncc.edu](mailto:rwhit101@uncc.edu)  
Or [open an issue](https://github.com/raw-lab/mercat2/issues)
