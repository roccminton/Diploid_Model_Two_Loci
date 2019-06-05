#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 14:04:24 2018

@author: larocca
"""

import numpy as np
import importlib
import copy

model_name = 'Families'


#import model specific functions such as birthrate, mutation operator and probabilities, etc.
spc = importlib.import_module('model_specifications.{}'.format(model_name))


def get_summerized_list(stats_by_time,trait1,trait2,index_low=0, index_high=None, step=1):
    """Returns a list with the number of the handed in key per time in the
    population between the index *index_high* and *index_low*. The handed in values
    *trait1* need to be or 'homo' or 'hetero' and the second *trait2* needs to be a list
    containing or '0', '1' or both representing the absence or presence of the recessive allele.
    
    """
    if index_high == None:
        index_high = len(stats_by_time)
    
    if trait1 not in ['homo','hetero']:
        raise ValueError('Handed in value in *get_summerized_list* must be or *hetero* or *homo*. {} is invalid.'.format(trait1))
    
    #get boolian for homogeneous or heterogeneous
    un_equal = (trait1=='homo')

    #generate empty list
    trait_per_time = []
    
    for count in stats_by_time[index_low:index_high:step]:
        #add a list entry
        trait_per_time.append(0)
        for trait in count['trait'].keys():
            #check if the trait is homogeneous/heterogeneous and if it carries the rezessive allele
            if (trait[0]==trait[1]) == un_equal and (trait[2]+trait[3]) in trait2:
                trait_per_time[-1] += count['trait'][trait]
        #check if in the current time interval there are any additions, if not set the value to None
        if trait_per_time[-1] == 0:
            trait_per_time[-1] = None
                
    return trait_per_time


    
def get_abs_quantities(stats_by_time,label_list,index_low=0, index_high=None, step=1):
    """Returns a list with the number of the handed in key per time in the
    population between the index *index_high* and *index_low*.
    
    The label must be one of (homo,n),(homo,n), for any integer n.    
    """
    for label in label_list:
        if label[0] not in ['homo','hetero'] or type(label[1]) != int or label[1] < 0 :
            raise KeyError('The handed in label {} in plotting_funcions.py is invalid'.format(label))
    
    if index_high == None:
        index_high = len(stats_by_time)

    trait_per_time = []
    
    for count in stats_by_time[index_low:index_high:step]:
        abs_quantity = 0
        for label in label_list:
            if label in count['pop_state'].keys():
                abs_quantity += count['pop_state'][label]
        if abs_quantity > 0:
            trait_per_time.append(abs_quantity)
        else:
            trait_per_time.append(None)
        
    return trait_per_time

def get_rel_quantities(stats_by_time,label_nominator,label_denominator,index_low=0, index_high=None, step=1):
    """Returns a list with the relative number of the handed in key per time in the
    population between the index *index_high* and *index_low* in all individuals with the key in the 
    denominator list. Note that both label need to be lists even if they content only one element.
    
    The label must be one of (homo,0),(homo,1),(hetero,0),(hetero,1).    
    """
    #Check if all handed in lables are valid
    for label in label_denominator + label_nominator:
        if label[0] not in ['homo','hetero'] or type(label[1]) != int or label[1] < 0 :
            raise KeyError('The handed in label {} in plotting_funcions.py is invalid'.format(label))
            
    #Check if the nominator lable is indeed a subgroup of the denominator lables
    #and afterwords delete it from the denominator list 
    #but save the orignial denominator in advance for later use
    old_denominator = copy.deepcopy(label_denominator)
    for label in label_nominator:
        if label in label_denominator:
            label_denominator.remove(label)
        else:
            raise KeyError('The handed in label {} in plotting_funcions.py is not part of the total quantity'.format(label))
    
    #check if nominator is a true subset of denominator
    if not label_denominator:
        raise KeyError('The handed in nominator label {} must be a true subset of the handed in denominator label {}'.format(label_nominator,old_denominator))

    
    if index_high == None:
        index_high = len(stats_by_time)

    trait_per_time = []
    
    for count in stats_by_time[index_low:index_high:step]:
        #Sum up nominator
        nominator = 0
        for label in label_nominator:
            if label in count['pop_state'].keys():
                nominator += count['pop_state'][label]
        #Sum up denominator starting with nominator and adding the rest
        denominator = nominator
        for label in label_denominator:
            if label in count['pop_state'].keys():
                denominator += count['pop_state'][label]
        if nominator > 0:
            trait_per_time.append(nominator/denominator)
        else:
            trait_per_time.append(None)
        
    return trait_per_time
    

def get_data(stats_by_time,data_type,index_low=0, index_high=None, step=1):
    """Returns a list with all data_type values in the count dictionary between
    index_low and index_high with given step size. Valid data_types are 'time',
    'trait','population_size','mutation_counter','mutation_load','number_false'
    and 'number_of_types'
    
    """
    valid_data_types = ['time','trait','population_size','mutation_counter','mutation_load','number_false','number_of_types']
    
    if data_type in valid_data_types == False:
        raise ValueError('The data type {} does not exist.'.format(data_type))
    
    if index_high == None:
        index_high = len(stats_by_time)
        
    data = []
    
    for count in stats_by_time[index_low:index_high:step]:
            data.append(count[data_type])
    
    return data

def get_data_over_pop_size(stats_by_time,data_type,index_low=0,index_high=None,step=1):
    """Returns a list with the quotient of all data_type_nom values over data_type_denom values in the count
    dictionary between index_low and index_high with given step size. Valid data_types are 'time',
    'trait','population_size','mutation_counter','mutation_load','number_false'
    and 'number_of_types'
    
    """
    valid_data_types = ['time','trait','population_size','mutation_counter','mutation_load','number_false','number_of_types']
    
    if data_type in valid_data_types == False:
        raise ValueError('The data type {} does not exist.'.format(data_type))
    
    if index_high == None:
        index_high = len(stats_by_time)
        
    data = []
    
    for count in stats_by_time[index_low:index_high:step]:
        data.append(count[data_type]/sum(count['pop_state'].values()))
    
    return data


def get_index_for_time(stats_by_time,t):
    """Returns the index of the list *stats_by_time* when first the time *t* is surpassed

    """
    if t==0:
        return 0

    for count in stats_by_time:
        if count['time'] >= t:
            return stats_by_time.index(count)

    raise ValueError('Handed in time in get_index_for_time is out of range.')

def set_time_index(t_0,t_max,stats_by_time):
    """Returns the index of t_0 and t_max or sets the indices to the lenght of stats_by_time
    or 0 if necessary.
    """
    if t_0 == 0:
        index_low = 0
    else:
        index_low = get_index_for_time(stats_by_time,t_0)
    
    if t_max ==  None:
        index_high = len(stats_by_time)
    else:
        index_high = get_index_for_time(stats_by_time,t_max)
        
    return index_low, index_high

def get_rel_entropy_of_family_sizes(stats_by_time,index_low=0, index_high=None, step=1):
    """Returns the entropy of the family sizes over time.
    
    """
    if index_high == None:
        index_high = len(stats_by_time)
        
    data = []
    
    for count in stats_by_time[index_low:index_high:step]:
        pop_size = sum(count['pop_state'].values())
        max_entropy = np.log2(pop_size)
        entropy = 0
        for size, freq in count['fam_distribution'].items():
            entropy += freq*(size/pop_size)*np.log(size/pop_size)
        data.append(-entropy/max_entropy)
        
    return data
            
    
            

