#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 13:04:37 2021

@author: nicolaskylilis
"""

def get_representatives(parent):
    """
    function integrates information from NCBI taxonomy resource, NCBI assembly resource, and NCBI genomes resource to return a table of species info for representative genomes accessions

    Parameters
    ----------
    parent : list of strings or ints
        name or tax_id of taxonomic branch.

    Returns
    -------
    df_clade_reference : pandas dataframe
        table of species with info of :['# assembly_accession', 'bioproject', 'biosample', 'wgs_master',
       'refseq_category', 'taxid', 'species_taxid', 'organism_name',
       'infraspecific_name', 'isolate', 'version_status', 'assembly_level',
       'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter',
       'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
       'excluded_from_refseq', 'relation_to_type_material',
       'asm_not_live_date', 'Replicons']

    """
    
    #### phylogenetic tree information: purpose here is to extract the NCBI tax ids of organisms under a tree branch
    
    # packages
    import ete3
    import pandas as pd
    
    # object: Provides a local transparent connector to the NCBI taxonomy database.
    ncbi = ete3.NCBITaxa()
    #ncbi.update_taxonomy_database()     # Updates the ncbi taxonomy database by downloading and parsing the latest taxdump.tar.gz file from the NCBI FTP site (via HTTP).
    
    df_clade_all = pd.DataFrame([])
   
    for entry in parent:
        # given a parent taxid or scientific species name, returns a list of all its descendants taxids. If intermediate_nodes is set to True, internal nodes will also be dumped.
        l_descendants = ncbi.get_descendant_taxa(entry, intermediate_nodes=False, rank_limit=None, collapse_subspecies= False, return_tree= False)
        dict_descendants = ncbi.get_taxid_translator(l_descendants)
        
        tree_descendants = ncbi.get_descendant_taxa(entry, intermediate_nodes=False, rank_limit=None, collapse_subspecies= False, return_tree= True)
        #print(tree_descendants.get_ascii(attributes=['sci_name', 'taxid']))
        
        df_clade = pd.DataFrame(dict_descendants.values(), index = dict_descendants.keys(), columns=["name"])
        df_clade.index.name = "tax_id"
        
        df_clade_all = df_clade_all.append(df_clade)
    
    del df_clade
    df_clade = df_clade_all
    
    ####  detect reference species from descendents: purpose here is to ectract representaive genomees and the genomes accession numbers
    
    # packages
    import pandas as pd
    import urllib.request
    import os
    
    if not os.path.exists("primary_data"):
        os.mkdir("primary_data")
        
    # assembly report refseq
    # Reference genomes: are genome assemblies that are annotated and updated by the assembly submitters and chosen by the RefSeq curatorial staff based on their quality and importance to the community as anchors for the analysis of other genomes in their taxonomic group.
    # Representative genomes: For species without a reference genome, one assembly is selected as representative among the live RefSeq assemblies 
    # more information @ https://www.ncbi.nlm.nih.gov/refseq/about/prokaryotes/#referencegenome
    link = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
    fpath = "primary_data/assembly_summary_refseq.txt"
    if not os.path.exists(fpath):
        print("\n***** MESSAGE *******\nFrom: pangenome_analysis_fcns.get_representatives()\nDownloading doc: assembly_summary_refseq.txt from NCBI FTP \n")
        urllib.request.urlretrieve(link, fpath, reporthook = show_progress)
    df_assembly_refseq= pd.read_csv(fpath, sep='\t', skiprows=(1), low_memory=False)
    df_assembly_refseq = df_assembly_refseq[df_assembly_refseq["assembly_level"] == 'Complete Genome']
    df_assembly_refseq = df_assembly_refseq[(df_assembly_refseq['refseq_category'] == 'representative genome') | (df_assembly_refseq['refseq_category'] == 'reference genome')]
    df_assembly_refseq.set_index('taxid', drop = False, inplace = True)
    
    # remove duplicated eentries for species tax id 
    dupl = df_assembly_refseq[df_assembly_refseq.duplicated(subset=['species_taxid'],  keep=False)]
    df_assembly_refseq = df_assembly_refseq[df_assembly_refseq["organism_name"] != "Escherichia coli O157:H7 str. Sakai"]
    
    # remove obligate endosymbionts
    df_assembly_refseq = df_assembly_refseq[df_assembly_refseq["organism_name"] != "Buchnera aphidicola str. Bp (Baizongia pistaciae)"]
   
    
    #### intersection of tree and reference genomes
    idx = (df_clade.index).intersection(df_assembly_refseq.index)
    df_clade_reference = df_assembly_refseq.loc[idx]
    
    # prokaryotes
    link = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/" + "prokaryotes.txt"
    fpath = "primary_data/prokaryotes.txt"
    if not os.path.exists(fpath):
        print("\n***** MESSAGE *******\nFrom: pangenome_analysis_fcns.get_representatives()\nDownloading doc: prokaryotes.txt from NCBI FTP. Large file - this might take a while\n")
        urllib.request.urlretrieve(link, fpath, reporthook = show_progress)
    df_prok= pd.read_csv(fpath, sep='\t', low_memory=False)
    df_prok = df_prok[df_prok['Status'] == "Complete Genome"]
    df_prok = df_prok[(df_prok['Reference'] == "REPR") | (df_prok['Reference'] == "REFR") ]
    df_prok.set_index('TaxID', drop = False, inplace = True)
    

    
    
    # adding replicon information 
    idx2 = (df_prok.index).intersection(df_clade_reference.index)
    df_clade_reference['Replicons'] = df_prok['Replicons'].loc[idx2]
    df_clade_reference.reset_index(drop=True, inplace=True)
    
    return df_clade_reference, dict_descendants, tree_descendants

def show_progress(a,b,c):
    '''''Callback function 
    @a:Downloaded data block 
    @b:Block size 
    @c:Size of the remote file 
    '''  
    per=100.0*a*b/c  
    if per>100:  
        per=100  
    print('%.2f%%' % per)
        
        
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
                print('UPDATE: FASTA sequences downloading:', r)
                os.system('ncbi-acc-download --out ' + r_path + ' --format fasta ' + str(r))
    
    print('UPDATE: FASTA sequences have been downloaded')
    
    

    
def annotate_with_prokka(d_species):
    """
    

    Parameters
    ----------
    d_species : python dict 
        d_species['organism name'].keys(['tax_id', 'replicons', 'species_dir_path', 'replicons_fpath']).

    Returns
    -------
    None.

    """
    
    #packages and modules
    import os
    
    # genome annotation with prokka
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
                cmd = "prokka" + " --outdir " + outdir_name + " --locustag " + str(replicon_fname) + " --prefix " + str(replicon_name) + " " + replicon_path
                os.system(cmd)
            else:
                print("Directory: " + outdir_name + " exists --> Skipping annotation")
    
















