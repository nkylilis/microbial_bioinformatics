#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 09:59:57 2021

@author: nicolaskylilis
"""

#%%

def retrieve_genome_annotations(d_org):
    """
    

    Returns
    -------
    d_org : python dictionary
        dict. of dictionaries with key: organism_name, values: dict of organism info ('tax_id', 'replicons', 'kegg_id' etc)

    """
    import os
    dir_name = "species_genomes" # directory that stores genome annotation files
    
    for key in d_org:
        
        #print(key)
        org_dir = os.path.join(dir_name,key)
        lst_dir2 = os.listdir(org_dir)
        lst_dir2 = [x for x in lst_dir2 if "_annotated" in x]
        #print(lst_dir2)
        
        fname = []
        for repl in lst_dir2:
            repl_dir = os.path.join(org_dir,repl)
            #print(repl_dir)
            fname += [os.path.join(repl_dir,y) for y in os.listdir(repl_dir) if ".faa" in y]
            
        d_org[key]["fpath_faa"] = fname
        
        
        # get aa sequences from .faa genome annotation files
        
        import pandas as pd 
        
        name = []
        aa_seq = []
        
        from Bio import SeqIO
        for f in fname:
            for record in SeqIO.parse(open(f),'fasta'):
            #print(record.id, record.seq)
                name += [record.id]
                aa_seq += [str(record.seq)]
            
        df = pd.DataFrame(data = aa_seq, index = name, columns = ["sequence"])
        d_org[key]["proteins"] = df
        
    return d_org


#%% Kegg orthology assignment

def kegg_orthology_assignemnt(d_org):
    """
    

    Parameters
    ----------
    d_org : TYPE
        DESCRIPTION.

    Returns
    -------
    d_org : TYPE
        DESCRIPTION.
    df_presence_absence_ko : TYPE
        DESCRIPTION.

    """
    """
    converts table of gene group presence absence naming from fasta sequence names regions to kegg orthology terms

    Parameters
    ----------
    d_org : python dictionary
        DESCRIPTION.
     
    Returns
    -------
    d_org : python dict
        addition of organism protein sequences from genome annotation files.
    df_presence_absence_ko : pandas dataframe
        DESCRIPTION. table of gene group presence absence with kegg orthology terms

    """
    print("\n *** Progress update: Starting kegg orthology annotation...\n")
    
    # retrieve table of proteome aa sequences from genome annotation files
    d_org = retrieve_genome_annotations(d_org)
    
    #### blastp genome annotations vs kegg gene sequences
    import kegg_db
    import pandas as pd
    
    for org in d_org:
        kegg_org = d_org[org]["kegg_id"]
        if kegg_org != "NA":
            #print(org)
            
            # get org protein aa seq
            df_genes = kegg_db.annotate_org_with_kegg_pathways(kegg_org, reannotate = False)
            if 'aaseq' not in list(df_genes.columns):
                df_genes = kegg_db.annotate_org_with_kegg_pathways(kegg_org, reannotate = True)
            
            # write a fasta file with aa_sequences from kegg db
            df_proteins = df_genes[df_genes["aaseq"].notna()]
            
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            from Bio import SeqIO
            
            records =[]
            for index,row in df_proteins.iterrows():
                #print(index,row)
                record = SeqRecord(Seq(row["aaseq"]), id=row["kegg_identifier"], description=row["gene_name"])
                records +=[record]
                del record
    
            SeqIO.write(records, "kegg_org_records.faa", "fasta")
            del records, df_proteins, df_genes
            
            # write a fasta file with aa_sequences from genome annotation file
            df_proteins = d_org[org]["proteins"]
            
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
            from Bio import SeqIO
            
            records =[]
            for index,row in df_proteins.iterrows():
                #print(index,row)
                record = SeqRecord(Seq(row["sequence"]), id=index)
                records +=[record]
                del record
    
            SeqIO.write(records, "annotation_records.faa", "fasta")
            del records
            
            # sequence local alignment
            print("\n***MESSAGE***\nFrom kegg_annotation.kegg_orthology_assignemnt()\nStarting sequence alignemnt for genome annotated CDSs of organism against the kegg gene db for this organism: " + org)
            
            """
            # makeblastdb
            application = "makeblastdb"
            dbtype = "-dbtype 'prot'"
            file_in = "-in kegg_org_records.faa"
            space = " "
            
            import os
            cmd = application + space + dbtype + space + file_in
            os.system(cmd)
            
            # blastp
            application = "blastp "
            query_file = "-query annotation_records.faa "
            db = "-db kegg_org_records.faa "
            output_file = "-out results.out "
            evalue_threshold = "-evalue 1e-6 "
            max_num_returned_hits = "-max_target_seqs 5 "
            outfmt = "-outfmt 6"
            
            import os
            cmd = application + query_file + db + output_file + evalue_threshold + max_num_returned_hits + outfmt
            #cmd  = "blastp -query annotation_records.faa -db kegg_org_records.faa -out results.out -outfmt 6 -max_target_seqs 1 -evalue 1e-6"
            os.system(cmd)
            print("Finished blastp for " + org)
            """
            
            # blastp
            application = "blastp "
            query_file = "-query annotation_records.faa "
            subject = "-subject kegg_org_records.faa "
            output_file = "-out results.out "
            evalue_threshold = "-evalue 1e-6 "
            max_num_returned_hits = "-max_target_seqs 1 "
            outfmt = "-outfmt 6"
            
            import os
            cmd = application + query_file + subject + output_file + evalue_threshold + max_num_returned_hits + outfmt
            #cmd  = "blastp -query annotation_records.faa -db kegg_org_records.faa -out results.out -outfmt 6 -max_target_seqs 1 -evalue 1e-6"
            os.system(cmd)
            print("Finished blastp for " + org)
            
            # gene mapper
            df_mapping_genes = pd.read_csv("results.out", sep='\t', header = None)
            names = df_mapping_genes.columns.to_list()
            new_names = ["query acc.ver","subject acc.ver", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"]
            name_dict = dict(zip(names,new_names))
            df_mapping_genes.rename(columns=name_dict, inplace = True)
            df_mapping_genes.set_index("query acc.ver", inplace = True, drop = True)
            
            # making sure only identical sequences are returned
            df_mapping_genes_100 = df_mapping_genes[df_mapping_genes["% identity"] == 100]
            df_mapping_genes_100 = df_mapping_genes_100[df_mapping_genes_100['% identity'] == 100 ]
            df_mapping_genes_100 = df_mapping_genes_100[df_mapping_genes_100['mismatches'] == 0 ]
            df_mapping_genes_100 = df_mapping_genes_100[df_mapping_genes_100['s. start'] == 1 ]
            df_mapping_genes_100 = df_mapping_genes_100[df_mapping_genes_100['alignment length'] == df_mapping_genes_100['s. end'] ]
            d_org[org]["gene_mapper"] = df_mapping_genes_100
            
            # remove blast files
            os.remove("annotation_records.faa")
            os.remove("kegg_org_records.faa")
            os.remove("results.out")
            #os.system("rm *.faa.*")
    
        
    #### kegg orthology assignment to gene presence/absence table
    
    import pandas as pd
    
    # read document into dataframe
    fname = "pangenome_analysis/output_with_alignment/gene_presence_absence.csv"
    df_presence_absence = pd.read_csv(fname, sep = ",")
    df_presence_absence.set_index("Gene", inplace = True, drop = True)
    del fname
    
    # assign kegg ontology (kegg pathways) to genes
    new_lst = []
    for org in df_presence_absence.columns.to_list():
        if "_combined" in org:
            new_lst += [org]
    df_presence_absence_copy = df_presence_absence[new_lst]
    
    df_presence_absence_ko = pd.DataFrame([], index = df_presence_absence_copy.index)
    for col in df_presence_absence_copy.columns.to_list():
        #print(col)
        org = col.replace("_combined","")
    
        # from pangenome analysis presence/absence file get genes for particular organism
        df_ = df_presence_absence_copy[[col]]
        # load the organism gene mapper(genome_annotaion:kegg_identifier)
        try:
            df_gene_mapper = d_org[org]["gene_mapper"]
        except:continue
        # kegg orthology for organism
        try:
            fname = d_org[org]["kegg_id"] + ".txt"
            fpath = os.path.join("/Volumes/Passport/SpeedyMicrobes/kegg_database/annotations",fname)
            df_ko = pd.read_csv(fpath, sep ="\t")
            df_ko.set_index("kegg_identifier", inplace = True)
        except: continue
    
        #
        ko_lst = []
        for index, rows in df_.iterrows():
            gene_id = rows[col]
            try:
                s_ = df_gene_mapper.loc[gene_id]
                ko = df_ko['kegg_pathways'].loc[s_["subject acc.ver"]]
            except: ko = "False"
            ko_lst +=[ko]
    
        
        df_presence_absence_ko[col] = ko_lst
    
    return d_org, df_presence_absence_ko


#%% Kegg pathways assignment Level A

def kegg_pathway_maps_analysis_level_a(d_org):
    
    #### kegg orthology assignemnt
    d_org, df_presence_absence_ko = kegg_orthology_assignemnt(d_org)
    
    #### kegg pathways assignement
    # load kegg pathways ontology file
    import kegg_db
    df_ontology = kegg_db.load_brite_pathways_ontology()
    
    
    # assign pathways to KO ids
    df_ontology = df_ontology[df_ontology["level_B"] != "Global and overview maps"]
    df_ontology.set_index("pathway_id", inplace = True, drop = True)
    path_ids = df_ontology.index
    
    import pandas as pd

    df_presence_absence_pathways_new = pd.DataFrame([]) # storage variable
    df_presence_absence_pathways_level_A = pd.DataFrame([]) # storage variable
    df_presence_absence_pathways_level_A_set = pd.DataFrame([]) # storage variable
    import ast
    for index, row in df_presence_absence_ko.iterrows():
        org_cols = df_presence_absence_ko.columns.to_list()
        path = []
        level_ = []
        for org in org_cols:
            lst_paths = ast.literal_eval(row[org])
            if lst_paths == False:
                path += ["False"]
                level_ += [["False"]]
                set_ = list(set([item for sublist in level_ for item in sublist]))
            else:
                idx = pd.Index(lst_paths).intersection(path_ids)
                result = idx.to_list()
                if len(result) == 0: 
                    path += ["False"]
                    level_ += [["False"]]
                    set_ = list(set([item for sublist in level_ for item in sublist]))
                else:  
                    path += [result]
                    level_result = []
                    for p in result:
                        level_result +=[df_ontology["level_A"].loc[p]]
                    level_ += [level_result]
                    set_ = list(set([item for sublist in level_ for item in sublist]))
                    
        
        s_ = pd.Series(data = path, index = org_cols, name = index)
        df_presence_absence_pathways_new = df_presence_absence_pathways_new.append(s_)
        
        s_a = pd.Series(data = level_, index = org_cols, name = index)
        df_presence_absence_pathways_level_A = df_presence_absence_pathways_level_A.append(s_a)
        
        # organisms collapse into a single group
        s_a_set = pd.Series(data = [set_], index = ["grouped"], name = index)
        df_presence_absence_pathways_level_A_set = df_presence_absence_pathways_level_A_set.append(s_a_set)
        
    
    
    #### preparing for  plotting
    
    # read document into dataframe
    fname = "pangenome_analysis/output_with_alignment/gene_presence_absence.csv"
    df_presence_absence = pd.read_csv(fname, sep = ",")
    df_presence_absence.set_index("Gene", inplace = True, drop = True)
    del fname
    
    # number of isolates
    df_results = df_presence_absence[['No. isolates']]
    # percentages
    perc = (df_presence_absence['No. isolates']/len(d_org.keys()))*100
    df_results["perc_occurance"] = perc
    # categorisation
    assign = []
    lst_perc = df_results["perc_occurance"].to_list()
    for percentage in lst_perc:
        if percentage > 94: assign += ["core"]
        elif (percentage <94 and percentage >14): assign += ["shell"]
        elif percentage < 15: assign += ["cloud"]
    df_results["pangenome"] = assign
    
    # adding pathways
    df_results = pd.concat([df_results, df_presence_absence_pathways_level_A_set], axis=1, join="outer")
    df_results = df_results.fillna(0)
    
    
    #### plotting
    all_ = []
    for region in ["core", "shell", "cloud"]:
        
        df = df_results[df_results["pangenome"] == region]
        df.drop(columns = ['No. isolates', 'perc_occurance', 'pangenome'], inplace = True)
    
        import itertools
        list2d = df.values.tolist()
        merged = list(itertools.chain(*list2d))
        list2d = merged
        merged = list(itertools.chain(*list2d))
    
        counts = dict()
        for i in merged:
            counts[i] = counts.get(i,0) + 1
        counts.pop("False")
        
        labels = list(counts.keys())
        values = [x/sum(list(counts.values()))*100 for x in list(counts.values())]
        
        df_ = pd.DataFrame(values, index = labels, columns =([region]))
        all_ += [df_]
    df_regions = pd.concat(all_, axis=1, join='inner')
    
    ax = df_regions.plot.barh(colormap='brg') # https://matplotlib.org/stable/tutorials/colors/colormaps.html
    ax.set_xlabel("% out of total map annotations per gene set (pangenome section)")
    ax.set_title("KEGG pathway maps - Level A")
    ax.set_axisbelow(True)
    ax.grid(color='gray', linestyle='dashed')
    
    #### save figure
    import os
    fname = "pangenome_kegg_annotation_pathways_level_A.png"
    dir_ = "results_plots"
    fpath = os.path.join(dir_,fname)
    import matplotlib.pyplot as plt
    plt.savefig(fpath, bbox_inches = 'tight') 
    
#%%  Kegg pathways assignment Level B

def kegg_pathway_maps_analysis_level_b(d_org, custom_list = []):
    """
    Description
    Creates bar plots figures of of level B pathway maps per pangenome section for specified Level A kegg categories. 
    Newly created figures are saved in dir: results_plots/
    *** Recomendation *** default custom list will take very long to compute. Use custom list based on results from the kegg pathway analysis for Level A maps

    Parameters
    ----------
    custom_list : list, optional
        DESCRIPTION. The default is []. 
                     Example custom_list = ['Cellular Processes','Human Diseases','Genetic Information Processing','Metabolism','Environmental Information Processing']
                     See all options by : import kegg_db
                                          df_ontology = kegg_db.load_brite_pathways_ontology(ontology_fpath = "primary_data/br08901.keg")
                                          print("Available Level A pathway maps options:\n",list(set(df_ontology["level_A"])))
                     
    Returns
    -------
    None.

    """
    
    #### kegg orthology assignemnt
    d_org, df_presence_absence_ko = kegg_orthology_assignemnt(d_org)
    
    
    #### pathway assignment
    # load kegg pathways ontology file
    import kegg_db
    df_ontology = kegg_db.load_brite_pathways_ontology(ontology_fpath = "primary_data/br08901.keg")
        
    if len(custom_list) == 0:
        level_A_lst = list(set(df_ontology["level_A"]))
    elif len(custom_list) != 0:
        level_A_lst = custom_list
        for x in custom_list:
            if x not in list(set(df_ontology["level_A"])):
                raise Exception(x + " is not a Level A kegg pathway map")
    
    
    import pandas as pd
    
    for level_A in level_A_lst:
        
        message = "*** UPDATE *** \nFrom kegg_annotation.kegg_pathway_maps_analysis_level_b(): \n - Annotating gene groups for Level A pathway maps subcategory: " + level_A + "\n..."
        print(message)
        
        # assign pathways to KO ids
        df_ontology_m = df_ontology[df_ontology["level_B"] != "Global and overview maps"]
        df_ontology_m = df_ontology_m[df_ontology_m["level_A"] == level_A]
        df_ontology_m.set_index("pathway_id", inplace = True, drop = True)
        path_ids = df_ontology_m.index
        
        
        df_presence_absence_ko_new = pd.DataFrame([]) # storage variable
        df_presence_absence_pathways_level_B = pd.DataFrame([]) # storage variable
        df_presence_absence_pathways_level_B_set = pd.DataFrame([]) # storage variable
        import ast
        for index, row in df_presence_absence_ko.iterrows():
            org_cols = df_presence_absence_ko.columns.to_list()
            path = []
            level_ = []
            for org in org_cols:
                lst_paths = ast.literal_eval(row[org])
                if lst_paths == False:
                    path += ["False"]
                    level_ += [["False"]]
                    set_ = list(set([item for sublist in level_ for item in sublist]))
                else:
                    idx = pd.Index(lst_paths).intersection(path_ids)
                    result = idx.to_list()
                    if len(result) == 0: 
                        path += ["False"]
                        level_ += [["False"]]
                        set_ = list(set([item for sublist in level_ for item in sublist]))
                    else:  
                        path += [result]
                        level_result = []
                        for p in result:
                            level_result +=[df_ontology_m["level_B"].loc[p]]
                        level_ += [level_result]
                        set_ = list(set([item for sublist in level_ for item in sublist]))
                        
            
            s_ = pd.Series(data = path, index = org_cols, name = index)
            df_presence_absence_ko_new = df_presence_absence_ko_new.append(s_)
            
            s_b = pd.Series(data = level_, index = org_cols, name = index)
            df_presence_absence_pathways_level_B = df_presence_absence_pathways_level_B.append(s_b)
            
            # organisms collapse into a single group
            s_b_set = pd.Series(data = [set_], index = ["grouped"], name = index)
            df_presence_absence_pathways_level_B_set = df_presence_absence_pathways_level_B_set.append(s_b_set)
        
        
            
        #### preparing for plotting
        
        # read document into dataframe
        fname = "pangenome_analysis/output_with_alignment/gene_presence_absence.csv"
        df_presence_absence = pd.read_csv(fname, sep = ",")
        df_presence_absence.set_index("Gene", inplace = True, drop = True)
        del fname
        
        # number of isolates
        df_results = df_presence_absence[['No. isolates']]
        # percentages
        perc = (df_presence_absence['No. isolates']/len(d_org.keys()))*100
        df_results["perc_occurance"] = perc
        # categorisation
        assign = []
        lst_perc = df_results["perc_occurance"].to_list()
        for percentage in lst_perc:
            if percentage > 94: assign += ["core"]
            elif (percentage <94 and percentage >14): assign += ["shell"]
            elif percentage < 15: assign += ["cloud"]
        df_results["pangenome"] = assign
        
        # adding pathways
        df_results = pd.concat([df_results, df_presence_absence_pathways_level_B_set], axis=1, join="outer")
        df_results = df_results.fillna(0)
    
        #### plotting    
            
        all_ = []
        for region in ["core", "shell", "cloud"]:
            
            df = df_results[df_results["pangenome"] == region]
            df.drop(columns = ['No. isolates', 'perc_occurance', 'pangenome'], inplace = True)
        
            import itertools
            list2d = df.values.tolist()
            merged = list(itertools.chain(*list2d))
            list2d = merged
            merged = list(itertools.chain(*list2d))
        
            counts = dict()
            for i in merged:
                counts[i] = counts.get(i,0) + 1
            counts.pop("False")
            
            labels = list(counts.keys())
            values = [x/sum(list(counts.values()))*100 for x in list(counts.values())]
            
            df_ = pd.DataFrame(values, index = labels, columns =([region]))
            all_ += [df_]
        df_regions = pd.concat(all_, axis=1, join='inner')
        
        ax = df_regions.plot.barh(colormap='brg') # https://matplotlib.org/stable/tutorials/colors/colormaps.html
        ax.set_xlabel("% out of total map annotations for level A map subcategory (per pangenome section)")
        ax.set_title("KEGG pathway maps - Level B:" + level_A)
        ax.set_axisbelow(True)
        ax.grid(color='gray', linestyle='dashed')
        
        #### save figure
        import os
        fname = "pangenome_kegg_annotation_pathways_level_B_" +level_A+".png"
        dir_ = "results_plots"
        fpath = os.path.join(dir_,fname)
        import matplotlib.pyplot as plt
        plt.savefig(fpath, bbox_inches = 'tight') 
    
    
    
        