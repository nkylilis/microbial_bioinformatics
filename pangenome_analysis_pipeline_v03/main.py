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
