#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 10:02:51 2021

@author: nicolaskylilis
"""

def get_doubling_data():
    """
    loads data file and produces table of species taxid and growth phenotype

    Returns
    -------
    df_phenotype : pandas dataframe
        index = species_tax_id; columns ['species', 'doubling_h'].

    """
    
    # packages and modules
    import pandas as pd
    
    # read database of phenotypic data
    df_phenotype = pd.read_csv("primary_data/Madin et al. (2020) - Bacterial and archaeal phenotypic trait data.csv", sep = ',')
    df_phenotype = df_phenotype[df_phenotype['doubling_h'].notna()]
    df_phenotype = df_phenotype[df_phenotype["growth_tmp"] >15]
    df_phenotype.set_index("species_tax_id",inplace=True)
    df_phenotype = df_phenotype[["species","doubling_h"]]
    
    return df_phenotype


