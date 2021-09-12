#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:33:41 2021

@author: nicolaskylilis
"""

#%% kegg db functions

def load_org_mapping():
    """
    loads a table to help map NCBI tax ids to kegg organism ids

    Returns
    -------
    df_organism_mapping : pandas dataframe
        returns a dataframe: index_col:ncbi_tax_id, columns:["kegg_tax_id", "name"]].

    """
    
    import pandas as pd
    df_organism_mapping = pd.read_csv("/Volumes/Passport/SpeedyMicrobes/kegg_database/kegg_organism_mapping.csv", sep = "\t", usecols =["ncbi_tax_id", "kegg_tax_id", "name"], index_col=("ncbi_tax_id"))
    
    return df_organism_mapping


    

def get_org_genes_from_kegg(org_kegg_identifier):
    """
    function retrieves all genes in kegg that have identifier thastarts with organism kegg id

    Parameters
    ----------
    org_kegg_identifier : str
        kegg db organism identifier ex. eco for E. coli MG1655.

    Returns
    -------
    df : pandas dataframe
        Index= 'kegg_identifier'(['gene_name'], dtype='object').

    """
    
    import Bio.KEGG.REST as REST
    
    response_str = REST.kegg_list(org_kegg_identifier).read()
    
    # read response string and make df
    line = []
    l_ = []
    for letter in response_str:
        if letter != "\n":
            line += letter
        elif letter == "\n":
            l_ +=["".join(line)]
            line=[]
    
    kegg_identifier = []
    gene_name = []
    for line in l_:
        l_items = line.split("\t")
        kegg_identifier += [l_items[0]]
        gene_name += [l_items[1]]
        
    
    
    import pandas as pd
    data = zip(kegg_identifier,gene_name)
    df = pd.DataFrame(data, columns=("kegg_identifier","gene_name"))
    df = df.set_index("kegg_identifier")
    
    return df


def annotate_org_with_kegg_pathways(org_kegg_identifier, reannotate = False): 
    """
    function annotates and save to file an org's genes on kegg_db with kegg pathways, kegg ko, NCBI_geneID, ncbi_proteID,
    
    Parameters
    ----------
    org_kegg_identifier : str
        kegg organism id.
    reannotate : boolean, optional
        reannotate the organsim genes to overwrite existing annotation file. The default is False.
        
    Returns
    -------
    df_genes: pandas dataframe.
        table of org genes labelled with ['kegg_identifier', 'gene_name', 'kegg_pathways', 'ncbi_geneID',
       'ncbi_protID', 'gene_location']
    
    """
    
    # check if file exists
    import os.path as path
    dir_path = "/Volumes/Passport/SpeedyMicrobes/kegg_database/annotations/"
    f_name = dir_path + org_kegg_identifier + ".txt"
    exists = path.exists(f_name)
    
    if (exists == True) & (reannotate == False): 
        print("\n*** MESSAGE***\nFrom: kegg_db.annotate_org_with_kegg_pathways()\n-Organism genes kegg pathways annotations loaded from file for: " + org_kegg_identifier + " \n")
        import pandas as pd
        df_genes = pd.read_csv(f_name, sep="\t")
    elif (exists == False) or (reannotate == True):
    
        import kegg_db
        df_genes = kegg_db.get_org_genes_from_kegg(org_kegg_identifier)
        #df_genes = df_genes.iloc[4500:]
        
        
        # find ko for genes in kegg database
        l_entries = []
        total_pathways = [] # storage variable
        total_geneID = [] # storage variable
        total_protID = [] # storage variable
        total_loc = [] # storage variable
        total_orth = []  # storage variable
        total_aaseq = []# storage variable
        count = 0
        total_gene_count = df_genes.shape[0]
        for KEGG_identifier, row in df_genes.iterrows():
            count +=1
            l_entries += [KEGG_identifier]
            
            #assign kos
            if len(l_entries) ==10:
                print(count, total_gene_count)
                #print("list of genes:",l_entries)
                [lst_geneID, lst_protID, lst_pathways, lst_loc, lst_orth, lst_aaseq] = kegg_db.get_gene_pathways(l_entries)
                
                # check that 10 parameters are returned
                if len(lst_geneID) + len(lst_protID) + len(lst_pathways) == (len(l_entries)*3):
                    pass
                else:
                    print("\n*** ERROR  ***\nFrom: kegg_db.annotate_org_with_kegg_pathways()\n - <10 response terms")
                    print("list of genes:",l_entries)
                    break
                
                #add results to storage variable
                total_pathways += lst_pathways
                total_geneID += lst_geneID 
                total_protID += lst_protID
                total_loc += lst_loc
                total_orth += lst_orth
                total_aaseq +=lst_aaseq
                
                #reset variable
                l_entries = []
                 
            elif (((df_genes.shape[0]-count)<10) & (df_genes.shape[0]==count)):
                [lst_geneID, lst_protID, lst_pathways, lst_loc, lst_orth, lst_aaseq] = kegg_db.get_gene_pathways(l_entries)
                
                #add results to storage variable
                total_pathways += lst_pathways
                total_geneID += lst_geneID 
                total_protID += lst_protID
                total_loc += lst_loc
                total_orth += lst_orth
                total_aaseq +=lst_aaseq
        
        
        df_genes["kegg_pathways"] = total_pathways
        df_genes["ncbi_geneID"] = total_geneID
        df_genes["ncbi_protID"] = total_protID
        df_genes["gene_location"] = total_loc
        df_genes["orthology"] = total_orth
        df_genes["aaseq"] = total_aaseq
        df_genes.reset_index(inplace = True)
        
        
        
        # save to file
        dir_path = "/Volumes/Passport/SpeedyMicrobes/kegg_database/annotations/"
        f_name = dir_path + org_kegg_identifier + ".txt"
        df_genes.to_csv(f_name, sep = "\t", index=False)
        
    return df_genes



def get_gene_pathways(gene_identifier):
    """
    get NCBI gene id, NCBI protein ID and kegg pathways for each kegg gene identifier by calling kegg's API

    Parameters
    ----------
    gene_identifier : list
        kegg db gene identifier lis ex.["eco:b0004", "eco:b0005", "eco:b0006", "eco:b0007", "eco:b0008", "eco:b0009", "eco:b00010"] .

    Returns
    -------
    lst_geneID : list
        list of NCBI protein_id lists for each identifier ex .['944771',
                                                                 '948295',
                                                                 '944751',
                                                                 '944750',
                                                                 '944753',
                                                                 '944754',
                                                                 '944756',
                                                                 '944758',
                                                                 '944757']
    lst_protID : list
        list of NCBI protein_id lists for each identifier ex .['NP_414552',
                                                                'YP_009518733',
                                                                'NP_414554',
                                                                'NP_414555',
                                                                'NP_414556',
                                                                'NP_414557',
                                                                'NP_414559',
                                                                'NP_414560',
                                                                'NP_414561']
    lst_pathways : list
        list of pathways lists for each identifier ex [['00260', '00750', '01100', '01110', '01120', '01230'],
                                                         False,
                                                         False,
                                                         False,
                                                         ['00030', '01100', '01110', '01120', '01200', '01230'],
                                                         ['00790', '01100', '01240', '04122']].
    lst_loc : list
        list of genes chromosome location ex ['1000:1050', '2000:3000']
        
    lst_aaseq: str

    """
    
    import Bio.KEGG.REST as REST
    response = REST.kegg_get(gene_identifier).read()
    #print(response)
    
    import kegg_db
    d_ = kegg_db.kegg_gene_parser(response)
    #print(list(d_['eco:b0004'].keys()))

    
    
    # get NCBI geneID 
    lst_geneID = []
    for key in d_.keys():
        gene_id = "False"
        try:
            l_ = (d_[key]["DBLINKS"])
            for x in l_:
                x = x.split(":")
                if x[0] =="NCBI-GeneID":
                    gene_id = x[1].strip()
                else: pass
        except: pass
        #print(key,gene_id)
        lst_geneID += [gene_id]
    #print(lst_geneID)
    
    
    # get NCBI protID 
    lst_protID = []
    for key in d_.keys():
        prot_id = "False"
        try:
            l_ = (d_[key]["DBLINKS"])
            for x in l_:
                x = x.split(":")
                if x[0] =="NCBI-ProteinID":
                    prot_id = x[1].strip()
                else: pass
        except:
            pass
        #print(key,prot_id)
        lst_protID += [prot_id]
    #print(lst_protID)
        
    
    # get pathways
    import re
    lst_pathways = []
    for key in d_.keys():
        try:
            l_ = d_[key]["PATHWAY"]
            l_new = []
            for x in l_:
                x = x.split(" ",1)[0]
                x = re.sub("[^0-9]", "", x)
                l_new += [x]
            lst_pathways += [l_new]
        except:
            lst_pathways += [False]
    
    # get gene location
    lst_loc = []
    for key in d_.keys():
        try:
            lst_loc += [d_[key]["POSITION"][0].replace("..",":")]
        except: lst_loc += [False]
    
    # get orthology
    lst_orth = []
    for key in d_.keys():
        try:
            lst_orth += [d_[key]["ORTHOLOGY"][0]]
        except:
            lst_orth += [False]
                
    
    # aa sequences
    lst_aaseq = []
    for key in d_.keys():
        try:
            lst_aaseq += [d_[key]["AASEQ"]]
        except:
            lst_aaseq += ["nan"]
        
    return lst_geneID, lst_protID, lst_pathways, lst_loc, lst_orth , lst_aaseq

def kegg_gene_parser(response):
    """
    
    Parameters
    ----------
    response : str
        str streeam from kegg's db API.
    
    Returns
    -------
    dict_gene:  dictionary
        keys: kegg entries from the str stream, values = ordered dict{}.
    
    """
    
    def next_item(odic, key):
        return list(odic)[list(odic.keys()).index(key) + 1]
    
    #packages and modules
    from collections import OrderedDict
    import random
    
    # storage variable
    dict_gene = OrderedDict() 
    
    
    l_response = response.split("///")[0:-1]
    # reult is text from response for individual gene entries
    for result in l_response:
        result = result+ "///\n"
        text_file = open("temp_ko.txt", "w")
        text_file.write(result)
        text_file.close()
        
        
        d_ = OrderedDict()
        
        # make dict keys for each row title of the particular kegg gene table with a value of the line that the text for that title begins
        count = 0
        with open("temp_ko.txt") as f_handle:
            for line in f_handle:
                #print(line)
                words = line.split(" ")
                #print(words[0])
                count +=1
                if words[0] != "":
                    #print(words[0])
                    d_[words[0]] = {}
                    d_[words[0]]["start"] = count
        
        # add the text for each row title [key] of the particular kegg gene table
        count = 0
        with open("temp_ko.txt") as f_handle:
            for line in f_handle:
                #print(line)
                count +=1
                words = line.split(" ")
                try: 
                    if words[0] in d_.keys():
                        if (count == d_[words[0]]["start"]):
                            line = line.replace(words[0],"")
                            line = line.strip()
                            d_[words[0]]["text"] = [line]
                            safe = words[0]
                    elif words[0] =="":
                        if ((count > d_[safe]["start"]) and (count < (d_[next_item(d_, safe)]["start"]))):
                            line = line.strip()
                            d_[safe]["text"] += [line]
                except:
                    pass
        # remove "///\n" from the dictionaryif exists   
        try:
            d_.pop("///\n")
        except:pass
        
        
        # copy the dict to make formating changes to keys(table titles) values
        dict_ = OrderedDict()       
        for key in d_.keys():
            dict_[key] = d_[key]["text"]
            
        #% Table title: ENTRY - fixed
        try:
            l_ = dict_["ENTRY"][0].split(" ")      
            dict_["ENTRY"] = [x for x in l_ if x != '']
        except: dict_["ENTRY"] = ["error_rand_" + str(random.randint(1,1000000000000000))]
            
        #% Table title: ORTHOLOGY - fixed
        try:
            lst = dict_["ORTHOLOGY"][0].split(" ",1)
            dict_["ORTHOLOGY"] = [x.strip() for x in lst]
        except: pass # dict_["ORTHOLOGY"] = ["error_rand_" + str(random.randint(1,1000000000000000))]
    
        #% Table title: DBLINKS - all good here
     
        #% Table title: PATHWAY - all good here  
    
        #% Table title: ORGANISM - fixed 
        try:
            lst = dict_["ORGANISM"][0].split(" ",1)
            dict_["ORGANISM"] = [x.strip() for x in lst]
        except: pass # dict_["ORGANISM"] = ["error_rand_" + str(random.randint(1,1000000000000000))]
        
        #% Table title: POSITION - fixed
        try:
            pos = dict_["POSITION"][0].replace("1:","")
            if "complement" in pos:
                pos = pos.replace("complement(","")
                pos = pos.replace(")","")
                l_pos = pos.split("..")
                pos = l_pos[1] + ".." + l_pos[0]
                dict_["POSITION"] = [pos]
            else:
                dict_["POSITION"] = [pos]
        except: pass #dict_["POSITION"] = ["error_rand_" + str(random.randint(1,1000000000000000))]
        
        #% Table title: AASEQ - fixed
        try:
            aa = dict_["AASEQ"][1:]
            aa = "".join(aa)
            dict_["AASEQ"] = aa
        except: dict_["AASEQ"] = None
        
        # storage variable
        try:
            dict_gene[dict_["ORGANISM"][0] + ":" + dict_["ENTRY"][0]] = dict_
        except:pass
        
    return dict_gene


def load_kegg_pathway_maps():
    """
    load to df all pathway maps of kegg database

    Returns
    -------
    df_maps : pandas dataframe
        df of index: kegg_map_code, columns = [map_name].

    """
    import pandas as pd
    
    f_name = "/Volumes/Passport/SpeedyMicrobes/kegg_database/kegg_pathways.txt"
    df_maps = pd.read_csv(f_name, sep="\t",  dtype=str)
    df_maps.set_index("kegg_map_code", inplace = True)
    
    return df_maps


def load_brite_pathways_ontology(ontology_fpath = ""):
    """
    Read kegg ontology for genes and proteins in brite hierarchy

    Returns
    -------
    df_ontology : pandas df
        Index(['level_A', 'level_b', 'level_C', 'pathway_id'], dtype='object')

    """
    
    if ontology_fpath == "":
        # file from https://www.biostars.org/p/70755/ or https://www.kegg.jp/kegg-bin/show_brite?br08901.keg
        fpath = "/Volumes/Passport/SpeedyMicrobes/kegg_database/br08901.keg"
    else:
        fpath = ontology_fpath
    
        
    text_copy =[]
    with open(fpath) as fhandle:
        text = fhandle.readlines()
            
        for line in text:
            if line[0] == "A":
                A = line[4:].replace("</b>","").strip()
            elif line[0] == "B":
                B = line[1:].strip()
            elif line[0] == "C":
                C = line[1:].strip()
                C1 = C.split("  ")[0]
                C2 = C.split("  ")[1]
                text_copy += [A + "\t" + B + "\t" + C2 + "\t" + C1]
            else: pass
    
    with open("temp.txt","w") as fhandle2:
        for element in text_copy:
            fhandle2.write(element + "\n")
    
    import pandas as pd
    df_ontology = pd.read_csv("temp.txt", sep = "\t", header = None, dtype=(str))
    mapper = dict(zip(df_ontology.columns.to_list(),["level_A","level_B","level_C","pathway_id"]))
    df_ontology.rename(columns = mapper, inplace = True)
    
    
    import os
    os.remove("temp.txt")
    
    return df_ontology


def get_organism_kegg_id(tax_id):
    """
    get organism kegg id

    Parameters
    ----------
    tax_id : int
        NCBI tax id.

    Raises
    ------
    Exception
        DESCRIPTION. More than 1 entries in db for tax id

    Returns
    -------
    kegg_tax_id : str
        organism kegg id from KEGG db.

    """
    
    df_kegg_org_mapping = load_org_mapping()
    
    try:
        df = df_kegg_org_mapping.loc[[tax_id]]
        print("*** MESSAGE ***\nFrom: kegg_db.get_organism_kegg_id()\n  -SUCCESS\n")
    except:
        import pandas as pd
        df = pd.DataFrame([["NA","NA"]], columns =["kegg_tax_id", "name"], index = [tax_id])
        print("*** WARNING ***\nFrom: kegg_db.get_organism_kegg_id()\n  -Organism NOT FOUND in KEGG database: NCBI_taxonomy_ID = ", tax_id, "\n")
        
    if not df.shape[0] ==1:
        raise Exception("\n*** WARNING ***\nFrom: kegg_db.get_organism_kegg_id()\n  -More than one matches in KEGG db for tax id: ", tax_id, "\n")
    else:
        kegg_tax_id = df["kegg_tax_id"].item()
        
    return kegg_tax_id


    

    
