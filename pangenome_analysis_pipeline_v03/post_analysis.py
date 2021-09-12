#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 11:12:13 2021

@author: nicolaskylilis
"""

#%% polt summary statistics

def plot_summary_statistics():
    fname = "pangenome_analysis/output_with_alignment/summary_statistics.txt"
    
    import pandas as pd
    df = pd.read_csv(fname, sep = "\t", header = None)
    
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn3
    
    core = df[2].iloc[0]
    soft_core = df[2].iloc[1]
    shell = df[2].iloc[2]
    cloud = df[2].iloc[3]
    
    # (size A group exclussive, size B group exclussive, intersection A & B, size C group exclussive,intersection C & A, intersection C & B, intersection C with A&B )
    labels=( 'core genes','shell genes','cloud genes',)
    #venn3(subsets = (cloud, 0, shell, 0, 0, 0, core + soft_core), set_labels=labels)
    venn3(subsets = (core + soft_core, shell, 0,cloud,0, 0, 0), set_labels=labels)
    plt.title("Pangenome gene pool distribution")
    
    import os
    dir_name ="results_plots"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    
    
    f_name = 'summ_stats_venn.png'
    f_path = os.path.join(dir_name,f_name)
    plt.savefig(f_path)
    plt.show()
    
    

#%% plot core convergence

def plot_converging_of_core():
    
    #### import data
    import pandas as pd
    f_path = "pangenome_analysis/output_with_alignment/number_of_conserved_genes.Rtab"
    df = pd.read_csv(f_path, sep = "\t", header = None)
    df.index.name = "iteration"
    df.columns = list(range(1,len(df.columns)+1))
    del f_path
    
    
    #### plot data
    import matplotlib.pyplot as plt
    import numpy as np
    x = list(df.columns)
    for t in list(df.index):
        y = list(df.iloc[t])
        plt.scatter(x,y,c = 'k', marker ='.')
    plt.xticks(np.arange(min(x), max(x)+1, 1.0))
    plt.xlabel("genomes added to set")
    plt.ylabel("number of conserved genes")
    #plt.show()
    del x, y, t
    
    
    #### make fitted model
    """
    from scipy.optimize import curve_fit
    import numpy as np
    
    # define model
    def func(x, k, t , w):
        # function from: Tettelin et al, 2005 PNAS (SI eq.1¨: Fc = κ􏰊c*exp[􏰈n􏰀􏰍/τc] + ω􏰒)
        # κc is the amplitude of the exponential decay
        # τc is the decay constant that measures the speed at which Fc (n) converges to its asymptotic value
        # Ω measures the size of the core genome for n → ∞
        y = (k * np.exp(-x/t)) + w
        return y
    
    # data
    df_stats = pd.DataFrame([])
    mean_data = []
    for col in df.columns:
        m = np.mean(list(df[col]))
        mean_data += [m]
    del col
        
    df_stats["genomes"] = df.columns
    df_stats["mean"] = mean_data
    del mean_data, m
    
    xData = df_stats["genomes"]
    yData = df_stats["mean"]


    # plot mean values
    #plt.scatter(xData, yData, color = 'r', marker = "x", label = 'mean values')
    #plt.legend()
    
    # perform regression
    initial_guess = [1,1,df_stats["mean"].iloc[-1]]
    popt, pcov = curve_fit(func, xData, yData, p0 = initial_guess)
    print("regression param: k,t,w",popt)
    
    # plot fitted model
    xModel = np.arange(1,20,0.01)
    plt.plot(xModel, func(xModel, *popt), 'r', label = 'fitted model')
    plt.legend()
    del xModel
    """
    # save plot
    import os
    dir_name ="results_plots"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    
    f_name = 'core_convergence_plot.png'
    f_path = os.path.join(dir_name,f_name)
    plt.savefig(f_path)
    
    #return popt, pcov


#%% plot new genes convergence

def plot_convergence_of_new_genes(regression_model = "heaps law"):
    
    regression_model = "heaps law"
    
    #### import data
    import pandas as pd
    f_path = "pangenome_analysis/output_with_alignment/number_of_new_genes.Rtab"
    df = pd.read_csv(f_path, sep = "\t", header = None)
    df.index.name = "iteration"
    df.columns = list(range(1,len(df.columns)+1))
    del f_path
    
    
    #### plot data
    import matplotlib.pyplot as plt
    import numpy as np
    x = list(df.columns)
    for t in list(df.index):
        y = list(df.iloc[t])
        plt.plot(x,y,c ='k', marker ='.', linestyle = "")
        #plt.scatter(np.log(x),np.log(y),c = 'k', marker ='.', linestyle = "")
        #plt.loglog(x,y,c = 'k', marker ='.', linestyle ="")
    plt.xticks(np.arange(min(x), max(x)+1, 1.0))
    #plt.xticks(np.arange(np.log(min(x)), np.log(max(x)+1, 1.0)))
    # plt.xlabel("genomes added to set")
    # plt.ylabel("number of new genes")
    #plt.show()
    del x, y, t
    
    
    #### fit model

    # define models
    def func_exp(x, k, t , theta):
        # function from: Tettelin et al, 2005 PNAS (SI eq.1¨: Fc = κ􏰊c*exp[􏰈n􏰀􏰍/τc] + ω􏰒)
        # κc is the amplitude of the exponential decay
        # τc is the decay constant that measures the speed at which Fc (n) converges to its asymptotic value
        # theta measures the number of specific genes for n → ∞ 
        y = (k * np.exp(-x/t)) + theta
        return y
    
    def func_power_law(ln_x, ln_k, a):
        # function from: Tettelin H, Riley D, Cattuto C, Medini D. Comparative genomics: the bacterial pan-genome. Curr Opin Microbiol. 2008 Oct;11(5):472-7. doi: 10.1016/j.mib.2008.09.006. PMID: 19086349.
        #import numpy as np
        #y = k * np.power(x,(-a))
        ln_y = ln_k - (a*ln_x)
        return ln_y
    
    
    # data
    df_stats = pd.DataFrame([])
    mean_data = []
    median_data = []
    for col in df.columns:
        mean = np.mean(list(df[col]))
        mean_data += [mean]
        median = np.median(list(df[col]))
        median_data += [median]

        
    df_stats["genomes"] = df.columns
    df_stats["mean"] = mean_data
    df_stats["median"] = median_data
    del mean_data, mean, median_data, median, col
    
    
    if regression_model == "heaps law":
        log_xData = (np.array(np.log(df_stats["genomes"]))).reshape((-1,1))
        log_yData = np.log(df_stats["median"])
        plt.plot(np.exp(log_xData), np.exp(log_yData), color = 'r', marker = "x", linestyle = "")
        
        # perform linear fit
        import numpy as np
        from sklearn.linear_model import LinearRegression
        model = LinearRegression().fit(log_xData, log_yData)
        r_sq = model.score(log_xData, log_yData)
        # print('coefficient of determination:', r_sq)
        # print('intercept:', model.intercept_)
        # print('slope:', model.coef_)
        
        # plot fitted model
        ln_k = model.intercept_
        k = np.exp(ln_k)
        a = -model.coef_[0]
        logxModel = np.arange(np.log(1),np.log(df_stats['genomes'].iloc[-1]) ,0.01) # xdata point
        xModel = np.exp(logxModel)
        yModel = np.exp(func_power_law(logxModel,ln_k, a))
        plt.plot(xModel, yModel, 'r')
        plt.ylabel("number of new genes added to pangenome")
        plt.xlabel("number of genomes analysed")
        param = "inferred param: k =" + str(round(k,2)) + ", alpha=" + str(round(a,2) )
        plt.legend(['median values','fitted model - log transformed y=k[x]^-a',param, "r_sq = " + str(round(r_sq,2)) ])
        
        
        
    elif  regression_model == "exponential":
        print("WARNING: Exponetial model fitting not implemented yet")
        
        # xData = df_stats["genomes"]
        # yData = df_stats["mean"]

        # # plot data
        # plt.scatter(xData, yData, color = 'r', marker = "x", label = 'mean values')
        # plt.legend()
        
        
        # # perform regression
        # from scipy.optimize import curve_fit
        # import numpy as np
        # initial_guess = [1,1,df_stats["mean"].iloc[-1]]
        # popt, pcov = curve_fit(func, xData, yData, p0 = initial_guess)
        # print("regression param: k,t,w",popt)
        
        # # plot fitted model
        # xModel = np.arange(1,20,0.01)
        # plt.plot(xModel, func(xModel, *popt), 'r', label = 'fitted model')
        # plt.legend()
        # del xModel

    
    # save analysis plots
    import os
    dir_name ="results_plots"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    
    f_name = 'new_genes_convergence_plot.png'
    f_path = os.path.join(dir_name,f_name)
    plt.savefig(f_path)
    


#%% plot pangenome convergence

def plot_converging_of_pangenome(regression_model = "heaps law"):
    
    
    #### import data
    import pandas as pd
    f_path = "pangenome_analysis/output_with_alignment/number_of_genes_in_pan_genome.Rtab"
    df = pd.read_csv(f_path, sep = "\t", header = None)
    df.index.name = "iteration"
    df.columns = list(range(1,len(df.columns)+1))
    del f_path
    
    
    #### plot data
    import matplotlib.pyplot as plt
    import numpy as np
    x = list(df.columns)
    for t in list(df.index):
        y = list(df.iloc[t])
        plt.scatter(x,y,c = 'k', marker ='.')
    plt.xticks(np.arange(min(x), max(x)+1, 1.0))
    #plt.xlabel("genomes added to set")
    #plt.ylabel("number of genes in pangenome")
    #plt.show()
    del t, x, y
    
    
    #### fit model

    # define model
    def func_exponential(x, theta, D, kc, w):
        # function from: Tettelin et al, 2005 PNAS (SI eq.4)
        # D is the average number of genes per sequenced genome
        y = D
        return y
    
    def func_power_law(x, k, gamma):
        # function from: Tettelin H, Riley D, Cattuto C, Medini D. Comparative genomics: the bacterial pan-genome. Curr Opin Microbiol. 2008 Oct;11(5):472-7. doi: 10.1016/j.mib.2008.09.006. PMID: 19086349.
        y = k *(np.power(x,gamma))
        return y
    
    # experimental data
    df_stats = pd.DataFrame([])
    mean_data = []
    median_data = []
    for col in df.columns:
        mean = np.mean(list(df[col]))
        mean_data += [mean]
        median = np.median(list(df[col]))
        median_data += [median]
    
    
    df_stats["genomes"] = df.columns
    df_stats["mean"] = mean_data
    df_stats["median"] = median_data
    df_stats["median_log"] = median_data
    del col, mean, mean_data, median, median_data
    
    if regression_model == "heaps law":
        
        xData = df_stats["genomes"]
        yData = df_stats["median"]
        plt.plot(xData, yData, color = 'r', marker = "x", linestyle = "")
        
        # perform curve fit
        from scipy.optimize import curve_fit
        bounds = ([0,0],[np.inf,1])
        popt, pcov = curve_fit(func_power_law, xData, yData, bounds=bounds)
        #print("optimised_parameters(k,gamma): ", popt)
        
        # plot fitted model
        xModel = np.arange(1,df_stats['genomes'].iloc[-1] ,0.01) # xdata point
        plt.plot(xModel, func_power_law(xModel, *popt), 'r')
        plt.ylabel("pangenome (number of genes)")
        plt.xlabel("number of genomes analysed")
        gamma = "inferred param gamma: k =" + str(round(popt[0],2)) + ", gamma=" + str(round(popt[1],2))
        plt.legend(['median values','fitted model - y=k[x]^gamma',gamma])
        #plt.show()
        
        
    elif regression_model == "exponential":
        # Warning: not tested thoroughly
        print("WARNING: exponential model fit not implemented yet")

        # xData = df_stats["genomes"]
        # yData = df_stats["mean"]
        # plt.plot(xData, yData, color = 'r', marker = "x", label = 'mean values & sd')
        
        # # perform curve fit
        # from scipy.optimize import curve_fit
        # popt, pcov = curve_fit(func_exponential, xData, yData)
        
        # # plot fitted model
        # xModel = np.arange(1,df_stats['genomes'].iloc[-1] ,0.01) # xdata point
        # plt.plot(xModel, func_exponential(xModel, *popt), 'r', label = 'fitted model, param: a=%5.3f, b=%5.3f' % tuple(popt))
        # plt.ylabel("pangenome(number of genes)")
        # plt.xlabel("number of genomes analysed")
        # plt.legend([])
        # plt.show()

        

    # save analysis plots
    import os
    dir_name ="results_plots"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    
    f_name = 'pangenome_convergence_plot.png'
    f_path = os.path.join(dir_name,f_name)
    plt.savefig(f_path)
    

#%% visualise gene absence/presence

def visualise_gene_presence_absence():

    import pandas as pd
    fpath = "pangenome_analysis/output_with_alignment/gene_presence_absence.Rtab"

    df= pd.read_csv(fpath, sep="\t")
    df.set_index("Gene", drop = True, inplace =True)
    df = df.transpose()


    import matplotlib.pyplot as plt
    import numpy as np
        
    labels = [x.replace("_combined","") for x in df.index.to_list()]
     
    fig = plt.figure()  # figure so we can add axis
    ax = fig.add_subplot(111)  # define axis, so we can modify
    ax.matshow(df, aspect = "auto")  # display the matrix
    #ax.set_xticks(np.arange(len(labels)))  # show them all!
    ax.set_yticks(np.arange(len(labels)))  # show them all!
    #ax.set_xticklabels(labels)  # set to be the abbv (vs useless #)
    ax.set_yticklabels(labels)  # set to be the abbv (vs useless #)
    
    
    # save analysis plots
    import os
    dir_name ="results_plots"
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
    
    f_name = 'presence_absence_matrix.png'
    f_path = os.path.join(dir_name,f_name)
    plt.savefig(f_path,bbox_inches='tight')
    #plt.show()


#%%

"""
import pandas as pd
fname = "pangenome_analysis/output_with_alignment/clustered_proteins"
df = pd.read_csv(fname, sep = ":", header = None)
df.set_index(0, drop =True, inplace = True)
df.rename(columns= {1 : "cluster"}, inplace = True)

lst = []
for i,cluster in df.iterrows():
    cluster = [x.strip() for x in cluster.item().split("\t")]
    lst += [cluster]
df["cluster"] = lst

del cluster, fname, i, lst

#%%
# build tree
from Bio import Phylo
import matplotlib.pyplot as plt

fname = "pangenome_analysis/output_with_alignment/accessory_binary_genes.fa.newick"
with open(fname) as f:
    text = f.read()
text = text.replace("_combined","")

file1 = open("temp_file.txt","w")
file1.write(text)
file1.close()

tree = Phylo.read("temp_file.txt", "newick")
tree.ladderize()  # Flip branches so deeper clades are displayed at top
fig = plt.figure()
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes,do_show=False)

# save analysis plots
import os
dir_name ="results_plots"
if not os.path.exists(dir_name):
    os.mkdir(dir_name)


f_name = 'treet.png'
f_path = os.path.join(dir_name,f_name)
fig.savefig(f_path)


"""


