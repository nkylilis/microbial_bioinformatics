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
    - Scoary: bacterial gwas
"""


#%% select group of organisms for pangenome analysis & store information in a dict = d_species

#### select organisms based on phenotypic data
import phenotype_fcns
df_phenotype = phenotype_fcns.get_doubling_data()


#### select clade of interest for pangenome analysis
# packages and modules
import pangenome_analysis_fcns
import re

clade = ["Vibrionales","Pasteurellales", "Enterobacterales", "Alteromonadales"]
df,dict_,tree_ = pangenome_analysis_fcns.get_representatives(clade)
#print(tree_.get_ascii(attributes=['sci_name', 'taxid']))

# Uncomment section to include an outgroup
"""
outgroup = 'Escherichia coli str. K-12 substr. MG1655'
if outgroup != None:
    df_outgroup, dict_outgroup,tree_outgroup = pangenome_analysis_fcns.get_representatives(outgroup)
    
    df = pd.concat([df,df_outgroup])
    del outgroup
"""

#### intersection of phenotypic data and phylogenetic tree data
df.set_index('species_taxid',inplace=True, drop=False)
idx_a = (df_phenotype.index).intersection(df.index)
df = df.loc[idx_a]
df["doubling_h"] = df_phenotype["doubling_h"].loc[idx_a]
df.set_index("taxid", drop=False, inplace=True)



#%%  organisms genomes annotation

# packages and modules
import pandas as pd
import os

#% build an organism dict with tax_id  and NCBI replicons refseq accession 

# storage variable
d_species = {}
for tx_id, sp in df.iterrows():
    
    name = sp['organism_name'].replace(' ','_')
    
    # dict key: organism name
    d_species[name] = {}
    
    # tax id
    d_species[name]['tax_id'] = tx_id
    
    # NCBI replicons refseq accession 
    l_repl = sp['Replicons'].split(';')
    repl = []
    for r in l_repl:
        try:
            m = re.search(':(.+?)/', r).group(1)
        except:
            m = re.search(':(.+?)/?$', r).group(1)
        repl += [m]
    d_species[name]['replicons'] = repl

del repl, r, name, sp, l_repl, m, tx_id


# get fasta sequences from NCBI nucleotides resource
import pangenome_analysis_fcns
pangenome_analysis_fcns.get_fasta_seq(d_species)

# constructing dict for paths for species directories and their replicons
cwd = os.getcwd() + '/species_genomes'
for root, dirs, files in os.walk(cwd):
   for name in dirs:
       try:
           d_species[name]['species_dir_path'] = [os.path.join(root, name)] 
       except:
           #print(name)
           pass

# creates a list of replicons paths and store in dict
for species in d_species:
    
    l_replicon_path = [] # storage var  
    for root, dirs, files in os.walk(d_species[species]['species_dir_path'][0]):
        for name in files:
            if "DS_Store" in name:
                pass
            elif ".fasta" in name:
                #print(name)
                l_replicon_path += [os.path.join(root, name)]
                
        d_species[species]['replicons_fpath'] = l_replicon_path
        
    del l_replicon_path

for species in d_species:
    print(species)


# genome annotation with prokka
import pangenome_analysis_fcns
pangenome_analysis_fcns.annotate_with_prokka(d_species)

del cwd, dirs, files, name, root, species

  
#%%  Clustering analysis

# packages & modules
import os


#### assemble single gff files for roary

# move gff files (from prokka annotation) to directory
for species in d_species:
    dir_path = d_species[species]['species_dir_path'][0] 
    if not os.path.exists(dir_path + '/gff_files'):
        os.mkdir(dir_path  + '/gff_files')
    os.system('cp ' + dir_path +'/*_annotated/*.gff ' + dir_path + '/gff_files')
    d_species[species]['gff_dir_path'] = dir_path + '/gff_files/'

    
for species in d_species:
    
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
                            if ((line.startswith('##')) & (not line == '##FASTA\n')):
                                pass
                            elif ((line.startswith('##')) & (line == '##FASTA\n')):
                                break
                            elif not line.startswith('##'):
                                f1.write(line)

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


#### homology groups search by roary

# move combined gff files to directory for roary
if not os.path.exists('pangenome_analysis'):
    os.makedirs('pangenome_analysis')

for species in d_species:
    dir_path_dirname = d_species[species]["species_dir_path"][0] + "/gff_files/"
    l_ = list(os.listdir(dir_path_dirname))
    for f in l_:
        if "combined" in f:
            #print("prinout: ", dir_path_dirname + f)
            os.system("cp " + dir_path_dirname + f +" pangenome_analysis/")



# Roary: Rapid large-scale prokaryote pan genome analysis
# software help:
# import os
# os.system("roary -h")
# Usage:   roary [options] *.gff


# run pangenome analysis
# running roary with multi-FASTA alignment of the core genes 

perc_identity = str("-i 50 ")       # minimum percentage identity for blastp [95]
perc_isolates = str("-cd 99 ")      # percentage of isolates a gene must be in to be core [99]
cmd = "roary -f pangenome_analysis/output_with_alignment " + perc_identity + perc_isolates + "-e -mafft pangenome_analysis/*.gff"
os.system(cmd)
del cmd

"""
question: why not evalue threshold for blastp based homology?
"""


#%% GWAS

#### set up environment 
import os
stout = os.environ["PATH"].split(":")
if "/opt/anaconda3/bin" not in stout:
    os.environ["PATH"] += ":/opt/anaconda3/bin"
print('environ["PATH"]:')  
print(os.environ["PATH"].split(":"))

#### input files

# input: clustering analysis results file
fpath_genes = "pangenome_analysis/output_with_alignment/gene_presence_absence.csv"
df_abspre = pd.read_csv(fpath_genes, sep =",")


# input: traits table file 
df_data = df
df_fast = df_data[df_data["doubling_h"] < 1]
df_slow = df_data[(df_data["doubling_h"] > 1) & (df_data["doubling_h"] < 5)]
df_data = pd.concat([df_slow, df_fast])

# produce traits file for scoary
lst_names = [x.replace(" ","_") +"_combined" for x in list(df_data["organism_name"])]
lst_names = [x.replace('(','') for x in lst_names]
lst_names = [x.replace(')','') for x in lst_names]

lst_trait_1 = []
for x in list(df_data["doubling_h"]):
    if x<2: 
        x =1
    else: x = 0
    lst_trait_1 += [x]

import pandas as pd
df_trait = pd.DataFrame(lst_trait_1, columns=(["fast_growth"]))
df_trait.index = lst_names

if not os.path.exists("primary_data"):
    os.mkdir("primary_data")
df_trait.to_csv("primary_data/traits.csv", sep=",")


fpath_traits = "primary_data/traits.csv"
df_trait = pd.read_csv(fpath_traits, sep =",")


#
if not os.path.exists("scoary_dir"):
    os.mkdir("scoary_dir")
cmd = "scoary -t "+fpath_traits+" -g"+fpath_genes + " -o scoary_dir"
os.system(cmd)


#%% process end notification
import email_nk
email_nk.send_email(script_name ="panGWAS")