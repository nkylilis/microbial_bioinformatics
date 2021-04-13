#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 15:48:54 2021

@author: nicolaskylilis
"""
"""
Description:
    - Script assigns NCBI taxonomy id to NCBI microbial genomes table (reference/representative)

Input files:
    - 'data/prokaryotes_ref_repr.csv' (NCBI reference and representative genomes table)

Output file:
    - 'data/bacteria_reference_genomes.csv' (input tabble with NCBI tax_id as the table index)
"""

#%%
# http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

import pandas as pd

df = pd.read_csv('data/prokaryotes_ref_repr.csv')
l_ref_names = list(df['#Organism Name'])
d_ref_txids = ncbi.get_name_translator(l_ref_names)

l_ids = []
for name in l_ref_names:
    l_ids += d_ref_txids[name]
df.index = l_ids
df.index.name = 'tax_id'

df.to_csv('data/bacteria_reference_genomes.csv', sep='\t', columns=None, header=True, index=True, index_label='tax_ids')
