MerCat2: python code for versatile k-mer counter and diversity estimator for database independent property analysis (DIPA) for multi-omic analysis
================================================

![GitHub Logo](doc/mercat_workflow.jpg)

  
Installing MerCat2: 
 - Available via Anaconda: Enable BioConda repo and run `conda install mercat2`  <br/>

Source
 - Download mercat_setup.py from bin  <br/>
 - Run pyhon mercat_setup.py to install all required dependencies  <br/>
 - Go to downloaded mercat2 folder and run commands <br/>

Usage <br/>

  -i I path-to-input-file <br/>
  -f F path-to-folder-containing-input-files <br/>
  -k K kmer length<br/>
  -n N no of cores [default = all]<br/>
  -c C minimum kmer count [default = 10]<br/>
  -pro run mercat on protein input file specified as .faa<br/>
  -p run prodigal on nucleotide assembled contigs. Must be one of ['.fa', '.fna', '.ffn', '.fasta','fastq']<br/>
  -h, --help show this help message<br/>

 
  
Citing Mercat
-------------
If you are publishing results obtained using MerCat2, please cite:



CONTACT
-------

Please send all queries to Mounika Ramapuram Naik &nbsp;&nbsp;      <a href="mailto:mramapur@uncc.edu?"><img src="doc/gmail.png" style="width:25px;height:25px"/>    </a> &nbsp; &nbsp;  <br /> 
Dr. Richard Allen White III &nbsp;&nbsp;   <a href="mailto:rwhit101@uncc.edu?"><img src="doc/gmail.png" style="width:25px;height:25px"/>      </a>
 <br />
Or [open an issue](https://github.com/raw-lab/cerberus/issues).

