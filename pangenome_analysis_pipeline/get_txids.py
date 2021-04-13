#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 14:11:28 2021

@author: nicolaskylilis
"""
# http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

#%%



name2taxid_1 = ncbi.get_name_translator(['Vibrionales'])
name2taxid_2 = ncbi.get_name_translator(['Vibrio natriegens'])

print(ncbi.get_rank(name2taxid_2['Vibrio natriegens']))
# printout: {691: 'species'}

lineage = ncbi.get_lineage(691)
print(lineage)
# printout: [1, 131567, 2, 1224, 1236, 135623, 641, 662, 717610, 691]

l_descendants = ncbi.get_descendant_taxa('Vibrio harveyi group', intermediate_nodes=False, rank_limit=None, collapse_subspecies=False, return_tree=False)
descendants = ncbi.get_descendant_taxa('Vibrio harveyi group', intermediate_nodes=False, rank_limit=None, collapse_subspecies=True, return_tree=True)
print(descendants.get_ascii(attributes=['sci_name', 'taxid', 'rank']))
l_names = ncbi.translate_to_names(descendants)
print(l_names)

ncbi.get_name_translator(['Vibrio sp. 99-8-1'])
print(ncbi.get_rank(2607602))

tree = ncbi.get_descendant_taxa('Vibrionales', collapse_subspecies=True, return_tree=True)
print(tree.get_ascii(attributes=['sci_name', 'taxid', 'rank']))


#%% get descendants taxids

# http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

parent = 'Vibrionaceae'
l_descendants = ncbi.get_descendant_taxa(parent, intermediate_nodes=False, rank_limit=None, collapse_subspecies= False, return_tree= False)
l_descendants = ncbi.get_descendant_taxa(parent, intermediate_nodes=False, rank_limit=None, collapse_subspecies= True, return_tree= False)

l_species = []
for sp in l_descendants:
    entry = ncbi.get_rank([sp])
    if entry[sp] != 'species':
        #print(entry)
        pass
    elif entry[sp] == 'species':
        l_species += [sp]
    else:
        print('What is this? :', entry)

