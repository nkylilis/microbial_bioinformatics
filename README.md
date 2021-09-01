# microbial_bioinformatics
Software tools for carrying out research in microbial bioinformatics

## Tools:
## > pangenome_analysis_pipeline
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

## > pangenome_analysis_pipeline_v03
* ### Description: 
  * Python script for carrying out pangenome analysis for microbial species
* ### Dependencies:
  * ETEtoolkit - ete-ncbiquery for NCBI taxonomy resource
  * prokka - microbial genome annotation (https://github.com/tseemann/prokka)
  * Roary - pangenome analysis pipeline (https://github.com/sanger-pathogens/Roary)
  * ete3 - A Python framework for the analysis and visualization of trees. (http://etetoolkit.org/)
* ### Usage:
  * executed by main.py file

## > panGWAS_v01
* ### Description: 
  * Python scripts for carrying out microbial Genome Wide Association Studies using pangenome analysis
* ### Dependencies:
  * ETEtoolkit - ete-ncbiquery for NCBI taxonomy resource
  * prokka - microbial genome annotation (https://github.com/tseemann/prokka)
  * Roary - pangenome analysis pipeline (https://github.com/sanger-pathogens/Roary)
  * Scoary - bacterial GWAS (https://github.com/AdmiralenOla/Scoary)
  * SciPy - Python-based ecosystem of open-source software for mathematics, science, and engineering.
  * ete3 - A Python framework for the analysis and visualization of trees. (http://etetoolkit.org/)
* ### Usage:
  * executed by main.py file
  * phenotypic data file should be placed in primary_data dir. Write a function that returns a pandas dataFrame of index = tax_id;
  * output: reulsts files in scoary_dir created upon scipt execution
