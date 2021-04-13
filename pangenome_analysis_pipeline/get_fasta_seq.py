#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 19:27:51 2021

@author: nicolaskylilis
"""

def get_fasta_seq(d_species):
    """
    

    Parameters
    ----------
    d_species : dictionary
        contains species names and NCBI nucleotide accession codes for their replicons

    Returns
    -------
    None.

    """
    
    import os
    
    for key in d_species.keys():
        sp_name = key
        sp_repl_lst = d_species[key]['replicons']
        for r in sp_repl_lst:
            #print(r)
            sp_dir = 'species_genomes/' + sp_name + '/'
            r_path = sp_dir + r +'.fasta'
            if not os.path.exists(sp_dir):
                os.makedirs(sp_dir)
                #print(sp_dir)
            if not os.path.exists(r_path):
                os.system('ncbi-acc-download --out ' + r_path + ' --format fasta ' + str(r))
    
    print('UPDATE: FASTA sequences have been downloaded')