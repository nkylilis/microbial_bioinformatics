#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 12:08:31 2021

@author: nicolaskylilis
"""
def get_represenatives(clade):
    """
    

    Parameters
    ----------
    clade : TYPE
        DESCRIPTION.

    Returns
    -------
    df_interest : TYPE
        DESCRIPTION.

    """
    
    
    # packages
    
    import pandas as pd
    from ete3 import NCBITaxa
    
    #%% load representative genomes table
    
    df_repr = pd.read_csv('data/bacteria_reference_genomes.csv', sep = '\t', index_col=('tax_ids'))
    
    
    #%% get descendants taxids
    
    # http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    
    parent = clade
    l_descendants = ncbi.get_descendant_taxa(parent, intermediate_nodes=False, rank_limit=None, collapse_subspecies= False, return_tree= False)
    #d_names = ncbi.get_taxid_translator(l_descendants)
    
    #%% detect reference species from descendents
    df_interest = pd.DataFrame([])
    for i in l_descendants:
        try:
            df_interest =df_interest.append(df_repr.loc[i])
        except:
            #print(ncbi.get_taxid_translator([i]))
            pass
        
    return df_interest

