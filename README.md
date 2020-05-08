# GeneFinder

DNA sequence from a sample strain of Salmonella is analyzed for any pathogenic (disease causing) genes. To do this, the program outputs snippets of DNA that are likely to be protein-coding genes—a process known as gene finding or gene prediction. 

# Usage
Make sure [python](https://www.python.org/downloads/) is installed on your device. 
As explained above, a sample data DNA sequence of Salmonella is chosen for this analysis. However, any DNA sequence can be ran by this program.

To start, clone the repo, and change files in /data folder to data or DNA sequence of your choice. 

Before running the program, be sure to load the data in bash:

```bash
$ from load import load_seq
$ dna = load_seq("./data/X73525.fa")
```
 
 Then run:
```bash
python genefinder.py
```

It might take longer than a few min for program to evaluate the results depending on the size of the data. 

Once the results(snippets of DNA) ouput, the snippets can be used in genetic search engine such as protein-BLAST to confirm whether or not the genes predicted by the program are close matches to known genes, and if so, whether their function is pathogenic. 

***protein-BLAST*** contains a database in which there are various strains of Salmonella bacterium and their genes have already been identified and logged, including strains that are known to cause diseases such as Typhoid fever (Salmonella enterica typhi).

