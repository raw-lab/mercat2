# tutorial here

Usage:
-----
 * -i I        path-to-input-file
 * -f F        path-to-folder-containing-input-files
 * -k K        kmer length
 * -n N        no of cores [default = all]
 * -c C        minimum kmer count [default = 10]
 * -pro        protein input file (.faa) 
 * -p          protein mode on nucleotide assembled contigs. Must be one of ['.fa', '.fna', '.ffn', '.fasta', '.fastq']
 * -h, --help  show this help message


By default mercat assumes that inputs provided is one of ['.fa', '.fna', '.ffn', '.fasta']

> Example: To compute all 3-mers, run `mercat2 -i test.fa -k 3 -n 8 -c 10 -p`          
 
 The above command:
* Runs prodigal on `test.fa`, then runs mercat on the resulting protein file.            
* Results are generally stored in input-file-name_{protein|nucleotide}.csv and input-file-name_{protein|nucleotide}_summary.csv  
   * `test_protein.csv` and `test_protein_summary.csv` in this example  
* `test_protein.csv` contains kmer frequency count, pI, Molecular Weight, and Hydrophobicity metrics for individual sequences.  
* `test_protein_summary.csv` contains kmer frequency count, pI, Molecular Weight, and Hydrophobicity metrics for all unique kmers across all sequences in `test.fa`
* `test_protein_diversity_metrics.txt` containing the alpha diversity metrics.
  
Other usage examples:
---------------------    
*  `mercat2 -i test.fna -k 3 -n 8 -c 10`  
   Run mercat on nucleotide input - one of ['.fa', '.fna', '.ffn', '.fasta','.fastq']
    
*   `mercat2 -i test.fna -k 3 -n 8 -c 10 -p`  
    Run prodigal on nucleotide input, generate a .faa protein file and run mercat on it
    
*   `mercat2 -i test.faa -k 3 -n 8 -c 10 -pro`  
    Run mercat on a protein input (.faa)

* All the above examples can also be used with  `-f input-folder` instead of `-i input-file` option
  -  Example:  `mercat2  -f /path/to/input-folder -k 3 -n 8 -c 10` --- Runs mercat on all inputs in the folder