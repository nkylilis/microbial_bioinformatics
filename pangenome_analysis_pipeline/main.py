#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 17:03:59 2021

@author: nicolaskylilis
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 13:34:59 2021

@author: nicolaskylilis
"""

"""
Pangenome analysis

- softwares:
    - ETEtoolkit: ete-ncbiquery
    - prokka: bacterial genome anotation
    - roary: pangenome analysis
    - FastTree: phylogentic tree
"""
#%% dependencies

# packages
import re
import os

# custom modules
from get_species_set import get_represenatives
from get_fasta_seq import get_fasta_seq


#%% select clade of interest for pangenome analysis
clade = "Vibrionaceae"
df = get_represenatives(clade)

#% find replicon accession for each species
d_species = {}
for tx_id, sp in df.iterrows():
    name = sp['#Organism Name'].replace(' ','_')
    d_species[name] = {}
    d_species[name]['tx_id'] = tx_id
    l_repl = sp['Replicons'].split(';')
    repl = []
    for r in l_repl:
        m = re.search(':(.+?)/', r).group(1)
        repl += [m]
    d_species[name]['replicons'] = repl

#%% get fasta sequences from NCBI nucleotides resource

get_fasta_seq(d_species)

#%%  Genome annotation

# constructing dict for paths for species directories and their replicons
cwd = os.getcwd() + '/species_genomes'
for root, dirs, files in os.walk(cwd):
   for name in dirs:
       try:
           d_species[name]['species_dir_path'] = [os.path.join(root, name)] 
       except:
           #print(name)
           pass

# creates a list of replicons paths and inserts in species dict
l_species = list(d_species.keys())
for species in l_species:
    l_replicon_path = []      
    for root, dirs, files in os.walk(d_species[species]['species_dir_path'][0]):
        for name in files:
            if "DS_Store" in name:
                pass
            elif ".fasta" in name:
                #print(name)
                l_replicon_path += [os.path.join(root, name)]
        d_species[species]['replicons_fpath'] = l_replicon_path
    del l_replicon_path


# Genome annotation
"""
# prokka - rapid bacterial genome annotation
# software help
# import os
# os.system("prokka --help")
# usage: prokka [options] <contigs.fasta>
"""

replicon_path  = []
for key in d_species:
    for replicon_path in d_species[key]['replicons_fpath']:
        replicon_fname = replicon_path.split("/")[-1]
        #print(key, replicon_fname)
        replicon_name = replicon_fname.replace(".fasta","")
        dir_name = replicon_path.replace(".fasta","")
        #os.system("prokka" + " --outdir " + dir_name + "_annotated"  + " --locustag " + str(replicon_fname) + " --prefix " + str(replicon_name) + " " + replicon_path)

        outdir_name = dir_name + "_annotated"
        if not os.path.exists(outdir_name):
            print("Directory: " + outdir_name + " doesnt exists --> Initiallising annotation pipeline")
            os.system("prokka" + " --outdir " + outdir_name + " --locustag " + str(replicon_fname) + " --prefix " + str(replicon_name) + " " + replicon_path)
        else:
            print("Directory: " + outdir_name + " exists --> Skipping annotation")
        #break
    #break
        
#%% assemble single gff files for bipartite genomes

# move gff files to directory
l_species = list(d_species.keys())
for species in l_species:
    dir_path = d_species[species]['species_dir_path'][0] 
    if not os.path.exists(dir_path + '/gff_files'):
        os.mkdir(dir_path  + '/gff_files')
    os.system('cp ' + dir_path +'/*_annotated/*.gff ' + dir_path + '/gff_files')
    d_species[species]['gff_dir_path'] = dir_path + '/gff_files/'

#%%      
for species in l_species:
    l_fnames = os.listdir(d_species[species]['gff_dir_path'])
    #print(l_fnames)
    # adding features
    for name in l_fnames:
        if '.gff' in name:
            if not "combined"in name:
                #print(name)
                with open(d_species[species]['gff_dir_path'] + name) as f:
                    with open(d_species[species]['gff_dir_path'] + species + '_combined.gff', "a") as f1:
                        for line in f:
                            if not line == '##FASTA\n':
                                f1.write(line)
                            elif line == '##FASTA\n':
                                break
    # adding fasta sequence
    for name in l_fnames:
        if '.gff' in name:
            if not "combined"in name:
                #print(name)
                with open(d_species[species]['gff_dir_path'] + name) as f:
                    with open(d_species[species]['gff_dir_path'] + species + '_combined.gff', "a") as f1:
                        logic = 0
                        count = 0
                        for line in f:
                            if line == '##FASTA\n':
                                logic = 1
                            if logic == 1:
                                if not line == '##FASTA\n':
                                    f1.write(line)


#%% pangenome analysis

# move combined gff files to directory
if not os.path.exists('pangenome_analysis'):
    os.makedirs('pangenome_analysis')

l_species = list(d_species.keys())
for species in l_species:
    dir_path_dirname = d_species[species]["species_dir_path"][0] + "/gff_files/"
    l_ = list(os.listdir(dir_path_dirname))
    for f in l_:
        if "combined" in f:
            #print("prinout: ", dir_path_dirname + f)
            os.system("cp " + dir_path_dirname + f +" pangenome_analysis/")
            
#%%   


"""
import os
if not os.path.exists('pang_analysis_gff_files'):
    os.mkdir("pang_analysis_gff_files")
    os.system("cp annotated_*/*.gff  pang_analysis_gff_files/")
"""

"""
# Roary: Rapid large-scale prokaryote pan genome analysis
# software help:
# import os
# os.system("roary -h")
# Usage:   roary [options] *.gff
"""

# run pangenome analysis
# running roary with multi-FASTA alignment of the core genes 

#blastp=95%
#os.system("roary -f genomes_ncbi_copy/pangenome_analysis/output_with_alignment -e -mafft genomes_ncbi_copy/pangenome_analysis/*.gff")
#blastp=50%
os.system("roary -f pangenome_analysis/output_with_alignment -i 50 -e -mafft pangenome_analysis/*.gff")

#%%
"""
# build phylogenetic tree with FastTree
import os
os.system("FastTree -nt -gtr output_with_alignment/core_gene_alignment.aln > tree.newick")
"""