# microbial_bioinformatics
Software tools for carrying out research in microbial bioinformatics

## Tools:
## 1) pangenome_analysis_pipeline
* ### Description: 
  * Python script for carrying out pangenome analysis for microbial species including genome annotation
* ### Dependencies:
  * ETEtoolkit - ete-ncbiquery for NCBI taxonomy resource
  * prokka - microbial genome annotation (https://github.com/tseemann/prokka)
  * Roary - pangenome analysis pipeline (https://github.com/sanger-pathogens/Roary)
  * FastTree - building phylogenetic trees ()
* ### Usage:
  * executed by main.py file
  * species fasta genome sequences to be downloaded in folders (folder name = species name) situated in the root directory
  * output: files in the pangenome_analysis directory created during script execution
