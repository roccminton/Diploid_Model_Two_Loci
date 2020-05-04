#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 11:53:28 2020

@author: larocca
"""
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
import os
import numpy as np
from sklearn.linear_model import LinearRegression
import plotting_functions as pf


#Remember to manually set the model name also in plotting_functions.py
model_name = 'Families'

#Absolute dir the script is in. Needs to be changed manually if the file is moveed
script_dir = r'/home/larocca/Dokumente/Simulationen/Diploid_Model_Two_Loci'

def load_data_and_return(mating_scheme,data_type):
    #set the relativ path
    rel_path = "out_analysis/data/averages_per_pop_size_nloci({}_{}).pickle".format(data_type,mating_scheme)
    
    # Load locations from pickle file
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "rb") as in_file:
        averages_dict = pickle.load(in_file)
    
    return averages_dict

def load_and_safe(mating_scheme,data_type):
    
    parameters_K_tmax = [
            (500,10000),
            (1000,15000),
            (2000,20000),
            (5000,25000),
            (10000,30000),
            (20000,35000),
            (50000,70000),
            (100000,100000)]
    
    parameters_n_loci = [50,300,500,1000,1500,2000,2500,3000,5000]

    
    total = len(parameters_K_tmax)*len(parameters_n_loci)
    run = 0
    
    averages = {}
    for K, t_max in parameters_K_tmax:
        for n in parameters_n_loci:
            #set the relativ path
            rel_path = "out_analysis/data/statistic_{}(K={},nloci={},mating={},n_runs=1".format(
                    model_name,
                    [K],
                    n,
                    mating_scheme,
                    )
            if mating_scheme == 'consang':
                rel_path += ",family_sizes=[100]).pickle"
            else:
                rel_path += ").pickle"
            # Load locations from pickle file
            abs_file_path = os.path.join(script_dir, rel_path)
            try:
                with open(abs_file_path, "rb") as in_file:
                    stats_by_time = pickle.load(in_file)
                    
                averages[(K,n)] = _get_average(pf.get_data_over_pop_size(stats_by_time,data_type),int(t_max/2),int(t_max))
            except FileNotFoundError:
                averages[(K,n)] = None
            
            #Print status to console
            run += 1
            print('Loading data: {}% ...please wait.'.format(int(100*(run/total))))
            
    print('Loading complete.')
            
    # Store list with locations after each round in a pickle file
    rel_path = "out_analysis/data/averages_per_pop_size_nloci({}_{})".format(data_type,mating_scheme)
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "wb") as out_file:
        pickle.dump(averages, out_file)

def _get_average(data,data_min,data_max):
    """Calculates the average of a list over an index range.
    
    """
    
    return sum(data[data_min:data_max])/len(data[data_min:data_max])
    

def plot(averages_dict,data_type,mating_scheme):
    
    if data_type == 'mutation_load':
        colors = 'OrRd_d'
    else:
        colors = 'GnBu_d'
    #setup the figure
    fig, ax = plt.subplots(figsize=(9.0,6.0))
    # Die obere und rechte Achse unsichtbar machen:
    ax.spines['top'].set_color('none')
    # Die untere Diagrammachse auf den Bezugspunkt '0' der y-Achse legen:
    ax.spines['bottom'].set_position(('data',0))
    #Titel setzten
    ax.set_xlabel("Number of Genes")
    
    if data_type == 'number_false':
        ax.set_ylabel('Equilibrium % of ill individual')
    elif data_type == 'mutation_load':
        ax.set_ylabel('Equilibrium Relative Mutation Load')
    else:
        raise ValueError('Add y-Axis lable for the new data_type *{}* before plotting.'.format(data_type))
        
    #set color map
    colors = sns.color_palette(colors, len(n_loci))
    ax.set_prop_cycle('color',colors)
    
    for size in pop_size:
        x = np.array(n_loci).reshape((-1,1))
        y = np.array([averages_dict[(size,n)] for n in n_loci])
        model = LinearRegression(n_jobs=-1).fit(x,y)
#        slopes[size] = model.coef_[0]
           
        ax.scatter(n_loci,y)
        ax.plot(n_loci,model.predict(x),label=size)
        
    #set the x-axis ticks
    #don't show the tick at pop_size 1000 because it is too close betwee 500-2000
    #and would be covered by the sorrounding ticks
    plt.xticks(ticks=n_loci,rotation=30)
    
    if data_type == 'number_false':
        ax.set_ylim(0,0.17)
        ax.set_yticks([0.025,0.05,0.075,0.1,0.125,0.15])
        ax.set_yticklabels(['2,5%','5%','7,5%','10%','12,5%','15%'])
    
    #show legend
    ax.legend(
            loc='upper left', 
            bbox_to_anchor=(0,1),
            fontsize=9,
            edgecolor='k',
            title = 'Equilibrium Population Size'
            )        
    
    path = "Reg_{}_vs_n_loci_{}.pdf".format(data_type,mating_scheme)
    #safe image 
    plt.savefig(os.path.join(r'/home/larocca/Dokumente/Simulationen/Diploid_Model_Two_Loci/out_analysis/figures',
                                     path))
    
    plt.show()

def plot_slopes(averages_dict,n_loci,pop_size):
    
    slopes = []
    for size in pop_size:
        x = np.array(n_loci).reshape((-1,1))
        y = np.array([averages_dict[(size,n)] for n in n_loci])
        model = LinearRegression(n_jobs=-1).fit(x,y)
        slopes.append(model.coef_[0])
    
    fig, ax = plt.subplots(figsize=(9.0,6.0))
    # Die obere und rechte Achse unsichtbar machen:
    ax.spines['top'].set_color('none')
    #Titel setzten
    ax.set_xlabel("Equilibrium Population Sizes")
    ax.set_ylabel("r")

    ax.plot(pop_size,slopes)
    
    plt.show()
    
def print_slopes(averages_dict,n_loci,pop_size):
    
    slopes = []
    for size in pop_size:
        x = np.array(n_loci).reshape((-1,1))
        y = np.array([averages_dict[(size,n)] for n in n_loci])
        model = LinearRegression(n_jobs=-1).fit(x,y)
        slopes.append(model.coef_[0])
    
    print(slopes)
    print(np.average(slopes))

#------------------------------------------------------------------------------

#set parameters
mating_scheme = 'random'
data_type = 'number_false'

if mating_scheme == 'consang':
    pop_size = [1000, 2000, 5000, 10000, 20000, 50000]
else:
    pop_size = [1000, 2000, 5000, 10000, 20000, 50000, 100000]
n_loci = [50, 300, 500, 1000, 1500, 2000, 2500, 3000, 5000]
    
averages = load_data_and_return(mating_scheme,data_type)
plot(averages,data_type,mating_scheme)