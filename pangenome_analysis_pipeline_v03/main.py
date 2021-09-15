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

"""
# make log file
f = open("my_log.log", "w")

#%% select group of organisms for pangenome analysis & store information in a dict = d_species

# packages and modules
import pangenome_analysis_fcns
import re

# select clade of interest for pangenome analysis
clade = "Vibrio"
df,dict_,tree_ = pangenome_analysis_fcns.get_representatives(clade)
print(tree_.get_ascii(attributes=['sci_name', 'taxid']))

# Uncomment section to include an outgroup
"""
outgroup = 'Escherichia coli str. K-12 substr. MG1655'
if outgroup != None:
    df_outgroup, dict_outgroup,tree_outgroup = pangenome_analysis_fcns.get_representatives(outgroup)
    
    df = pd.concat([df,df_outgroup])
    del outgroup
"""

#% build an organism dict with tax_id  and NCBI replicons refseq accession 

# storage variable
d_species = {}
for i, sp in df.iterrows():
    
    name = sp['organism_name'].replace(' ','_')
    
    # dict key: organism name
    d_species[name] = {}
    
    # tax id
    d_species[name]['tax_id'] = sp["taxid"]
    
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

del repl, r, name, sp, l_repl, m, i


#%%  organisms genomes annotation

#packages & modules
import os

# get fasta sequences from NCBI nucleotides resource
import pangenome_analysis_fcns
pangenome_analysis_fcns.get_fasta_seq(d_species)


# constructing dict for paths for species directories and their replicons
if not os.path.exists("species_genomes"):
    os.mkdir("species_genomes")
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

  
#%% assemble single gff files for roary

for species in d_species.keys():
    
    gff3 = [] # combined storage variable
    gff = []
    fasta = []
    
    for repl in d_species[species]["replicons"]:
            
        f_path_fasta = d_species[species]['species_dir_path'][0] + "/" + repl + ".fasta"
        with open(f_path_fasta) as fhandle2:
            fasta_txt = fhandle2.readlines()
            line_0 = fasta_txt[0].split(" ")[0] + "\n"
            fasta_txt = [line_0] + fasta_txt[1:-1]
        fasta += fasta_txt
            
        
        f_path_gff = d_species[species]['species_dir_path'][0] + "/" + repl + "_annotated" + "/" + repl + ".gff"
        with open(f_path_gff) as fhandle:
            gff_txt = fhandle.readlines()[1:]
            gff_txt_new =[]
            for line in gff_txt:
                if line[0:7] != "##FASTA":
                    gff_txt_new +=[line]
                else:
                    break
        gff += gff_txt_new

    gff3 = ["##gff-version 3\n"] + gff + fasta
    
    
    if not os.path.exists(d_species[species]['species_dir_path'][0] + "/gff_files"):
        os.mkdir(d_species[species]['species_dir_path'][0] + "/gff_files")
    fname = d_species[species]['species_dir_path'][0] + "/gff_files/" + species + "_combined.gff" 
    with open(fname, "w") as fhandle3:
        for item in gff3:
            fhandle3.write(item)

#%% pangenome analysis

# packages & modules
import os

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
#blastp=50%rroary
i = str(50) #minimum percentage identity for blastp

cmd = "roary -f pangenome_analysis/output_with_alignment -i " + i + " -e -mafft pangenome_analysis/*.gff"
os.system(cmd)
del cmd

#%% pangenome structure analysis

import post_analysis

# plotting results
# graphs are saved in dir: results_plots/
post_analysis.plot_summary_statistics()
post_analysis.plot_converging_of_core()
post_analysis.plot_convergence_of_new_genes()
post_analysis.plot_converging_of_pangenome()

#%% kegg pathways analysis of pangenome

import kegg_db

# get organism ids
for org in d_species:
    #print(key)
    tax_id = d_species[org]["tax_id"]
    print("\n Searching KEGG db for: " + org)
    d_species[org]["kegg_id"] = kegg_db.get_organism_kegg_id(tax_id)


# LEVEL A analysis
import kegg_annotation
kegg_annotation.kegg_pathway_maps_analysis_level_a(d_species)


# LEVEL B analysis
import kegg_annotation
kegg_annotation.kegg_pathway_maps_analysis_level_b(d_species, custom_list = ['Cellular Processes','Human Diseases','Genetic Information Processing','Metabolism','Environmental Information Processing'])
