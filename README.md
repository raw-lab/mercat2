MerCat2: python code for versatile k-mer counter and diversity estimator for database independent property analysis (DIPA) for multi-omic analysis
================================================

![GitHub Logo](doc/mercat_workflow.jpg)

  
Installing MerCat2: 
 - Available via Anaconda: Enable BioConda repo and run `conda install mercat2`  <br/>

**Source** **Installer**
 - Download mercat_setup.py from bin  <br/>
 - Run pyhon mercat_setup.py to install all required dependencies  <br/>
 - Go to downloaded mercat2 folder and run commands <br/>

**Usage** <br/>
 - -i I path-to-input-file <br/>
 - -f F path-to-folder-containing-input-files <br/>
 - -k K kmer length<br/>
 - -n N no of cores [default = all]<br/>
 - -c C minimum kmer count [default = 10]<br/>
 - -pro run mercat on protein input file specified as .faa<br/>
 - -p run prodigal on nucleotide assembled contigs. Must be one of ['.fa', '.fna', '.ffn', '.fasta','fastq']<br/>
 - -h, --help show this help message<br/>

By default mercat assumes that inputs provided is one of ['.fa', '.fna', '.ffn', '.fasta'] <br/>

       Example: To compute all 3-mers, run python mercat2.py -i RW2.fna -k 3 -n 8 -c 10 -p
 
The above command:

- Runs prodigal on test.fna, then runs mercat on the resulting protein file.<br/>
- Results are generally stored in input-file-name_{protein|nucleotide}.csv and input-file-name_{protein|nucleotide}_summary.csv
       - RW2_protein.csv and RW2_protein_summary.csv in this example
- RW2_protein.csv contains kmer frequency count, pI, Molecular Weight, and Hydrophobicity metrics for individual sequences.
- RW2_protein_summary.csv contains kmer frequency count, pI, Molecular Weight, and Hydrophobicity metrics for all unique kmers across all sequences in test.fna
- RW2_protein_diversity_metrics.txt containing the alpha diversity metrics.

**Other** **usage** **examples**:

- python mercat2.py -i test.faa -k 3 -n 8 -c 10 -pro</br>
     Run mercat on a protein input (.faa)
- python mercat2.py -i RW2.fna -k 3 -n 8 -c 10 -p </br>
      Run prodigal on nucleotide input, generate a .faa protein file and run mercat on it
- python mercat2.py -i RW2.fna -k 3 -n 8 -c 10 </br>
      Run mercat on nucleotide input - one of ['.fa', '.fna', '.ffn', '.fasta','fastq']

All the above examples can also be used with -f input-folder instead of -i input-file option</br>

      Example: python mercat2.py -f /path/to/input-folder -k 3 -n 8 -c 10 --- Runs mercat on all inputs in the folder

For folder input , PCA analysis for all the samples is plotted in PCA_plot.html.

PCA analysis result for 5 sample files in data:
doc/PCA.png
  
Citing Mercat
-------------
If you are publishing results obtained using MerCat2, please cite:



CONTACT
-------

Please send all queries to Mounika Ramapuram Naik &nbsp;&nbsp;      <a href="mailto:mramapur@uncc.edu?"><img src="doc/gmail.png" style="width:25px;height:25px"/>    </a> &nbsp; &nbsp;  <br /> 
Dr. Richard Allen White III &nbsp;&nbsp;   <a href="mailto:rwhit101@uncc.edu?"><img src="doc/gmail.png" style="width:25px;height:25px"/>      </a>
 <br />
Or [open an issue](https://github.com/raw-lab/cerberus/issues).

