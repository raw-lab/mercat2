MerCat2: python code for versatile k-mer counter and diversity estimator for database independent property analysis (DIPA) for multi-omic analysis
================================================

![GitHub Logo](doc/mercat_workflow.jpg)

  
**Installing MerCat2:** 
 - Available via Anaconda: Enable BioConda repo and run `conda install mercat2`  <br/>

**Source Installer**
 - Download mercat_setup.py from bin  <br/>
 - Run python mercat_setup.py to install all required dependencies  <br/>
 - Go to downloaded mercat2 folder and run commands <br/>

**Usage** <br/>
 - -i I path-to-input-file <br/>
 - -f F path-to-folder-containing-input-files <br/>
 - -k K kmer length<br/>
 - -n N no of cores [default = all]<br/>
 - -c C minimum kmer count [default = 10]<br/>
 - -pro run mercat on protein input file specified as .faa<br/>
 - -p run prodigal on nucleotide assembled contigs. 
    - Must be one of ['.fa', '.fna', '.ffn', '.fasta','fastq']<br/>
 - -h, --help show this help message<br/>

By default mercat assumes that inputs provided is nucleotide mode of ['.fa', '.fna', '.ffn', '.fasta'] <br/>

**Usage examples**:

****Run mercat2 on a protein mode (protein fasta - '.faa')****</br>
`python mercat2.py -i test.faa -k 3 -n 8 -c 10 -pro`</br>

****Run mercat2 on a nucleotide mode (nucleotide fasta - '.fa', '.fna', '.ffn', '.fasta')****</br>
`python mercat2.py -i RW2.fna -k 3 -n 8 -c 10 -p` </br>

****Run mercat2 on a nucleotide mode raw data (nucleotide fastq - '.fastq')****</br>
`python mercat2.py -i RW2.fastq -k 3 -n 8 -c 10 -p` </br>

****Run on many samples within a folder****</br>
`python mercat2.py -f /path/to/input-folder -k 3 -n 8 -c 10`</br>

**Outputs**
- Results are stored in input-file-name_{protein|nucleotide}.csv and input-file-name_{protein|nucleotide}_summary.csv </br>
   -  file-name_protein.csv and file-name_protein_summary.csv (for example) </br>
- file-name_protein.csv (If run in protein mode)</br>
   -  Contains kmer frequency count, pI, Molecular Weight, and Hydrophobicity metrics for individual sequences.</br>
- file-name_protein_summary.csv (If run in nucleotide mode)</br>
   -  Contains kmer frequency count, pI, Molecular Weight, and Hydrophobicity metrics for individual sequences.</br>
- file-name_diversity_metrics.txt </br>
   -  Contains the alpha diversity metrics.</br>
- If 5 samples or more it will generate a PCA, PCA analysis for all the samples is plotted in PCA_plot.html.</br>

PCA analysis result for 5 sample files in data:

 ![GitHub Logo](doc/PCA.png)
  
Citing Mercat
-------------
If you are publishing results obtained using MerCat2, please cite:


CONTACT
-------

Please send all queries to Mounika Ramapuram Naik &nbsp;&nbsp;      <a href="mailto:mramapur@uncc.edu?"><img src="doc/gmail.png" style="width:25px;height:25px"/>    </a> &nbsp; &nbsp;  <br /> 
Dr. Richard Allen White III &nbsp;&nbsp;   <a href="mailto:rwhit101@uncc.edu?"><img src="doc/gmail.png" style="width:25px;height:25px"/>      </a>
 <br />
Or [open an issue](https://github.com/raw-lab/mercat2/issues).

