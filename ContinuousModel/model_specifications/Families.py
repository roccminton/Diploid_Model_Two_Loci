"""The model specifications are due to an example from Dieckmann and Doebeli
    (1999, :cite:'Dieckmann1999'). Other specifications can be added by creating
    a corresponding .json and .py file in the IN_MODEL_SPECS folder.

"""

import numpy as np
import re

    
def trait_mutation_probability(model):
    """Calculates the markermutation in the given setup.

    """
    if type(model['inv_mutation_rate']) == str: 
        word = model['inv_mutation_rate']
        regexp = re.compile(r'K\^\(\d|\d\)')
        
        if word == 'per bp':
            return model['mutation_rate']
        elif word == 'K':
            return (model['K'])**(-1)
        elif word == 'Klog(K)':
            return (model['K']*np.log(model['K']))**(-1) 
        elif word == 'Klog^2(K)':
            return (model['K']*np.log(model['K'])**2)**(-1)
        elif regexp.search(word):
            return model['K']**(-float(word[3])/float(word[5]))
        else:
            raise ValueError(
                    "Unknown inverse mutation rate {}. Add the mutation"
                    "rate to the specification file and try again.".format(model['inv_mutation_rate'])
                    )
    elif type(model['inv_mutation_rate']) == float or type(model['inv_mutation_rate']) == int:
        return 1/model['inv_mutation_rate']
    else:
        raise TypeError("Unknow type for inverse mutation rate. {} is as {} not valid.".format(model['inv_mutation_rate'],type(model['inv_mutation_rate'])))
    
def trait_mutation(parent_trait,model):
    """Simulates a new mutant trait given the parent trait and returns the corresponding
    trait value.
    
    """
    #here mutation always gives the unfit rezessive trait
    if model['gene_distribution'] == 'unif':
        return 1
    else:
        return parent_trait + 1

def birth_rate(trait,model):
    """Calculates the birthrate for an individuum with given trait.

    """
    if trait[2]:
        return 1/54.56
    else:
        return 0

def death_rate(trait,model):
    """Calculates the birthrate for an individuum with given trait.

    """

    return 1/72.7

def rep_comp(x,y,model):
    """Sets the reproductive compatibility of two individuals with genotypes x and y.
    Here we want individuals with equal types to reproduce more often than these with
    different types.
    
    """
    #DYNASTY
    if model['mating_scheme'][0] == 'dynasty':
        for family in x[0:2]:
            #at least one match (aa-aa, aa-ab, ab-ab, ab-ac)
            if family in y[0:2]:
                return 1
        #no match
        return model['mating_scheme'][1]
   
    #WESTERN
    elif model['mating_scheme'][0] == 'western':
        for family in x[0:2]:
            #at least one match (aa-aa, aa-ab, ab-ab, ab-ac)
            if family in y[0:2]:
                return model['mating_scheme'][1]
        #no match
        return 1
    
    #RANDOM
    elif model['mating_scheme'][0] == 'random':
        return 1
    
    #CONSANGUIN
    elif model['mating_scheme'][0] == 'consang':
        #perfect match
        if x[0:2] == y[0:2]:
            return model['mating_scheme'][1]
        #one homogeneous and one half
        elif (x[0]==x[1] or y[0]==y[1]) and (x[0]==y[0] or x[0]==y[1] or x[1]==y[0]):
            return model['mating_scheme'][2]
        #at least one match
        elif x[0]==y[0] or x[0]==y[1] or x[1]==y[0] or x[1]==y[1]:
            return model['mating_scheme'][3]
        #no match at all
        else:
            return model['mating_scheme'][4]
        
    #%IN_OUT_FAMILY
    elif model['mating_scheme'][0] == '%in_out_family':
        if x[0:2] == y[0:2]:
            return model['mating_scheme'][1]/model['number_of_families_list'][model['index_change']]
        else:
            return model['mating_scheme'][2]/(model['K_list'][model['index_change']]-model['number_of_families_list'][model['index_change']])

    else:
        raise KeyError('*{}* is no valid mating scheme in the *.json file'.format(model['mating_scheme'][0]))
    
    
def competition(x,y,model):
    """ Calculates the competition pressure executed from an individual with 
    trait x to an individual with trait y. """
    
    #uniform competition pressure to get pop equilibrium at K
    return birth_rate(y,model)-death_rate(y,model)

def pop_equilibrium(trait,model,trait2=None):
    """Calculates the equilibrium population size for a population monomorphic
    with given trait in respect to birth- and deathrate plus competition.

    """
    if trait2==None:
        return (birth_rate(trait,model) - death_rate(trait,model))/(competition(trait,trait,model))
    else:
        r_x = birth_rate(trait,model) - death_rate(trait,model)
        r_y = birth_rate(trait2,model) - death_rate(trait2,model)
        c_xx = competition(trait,trait,model)
        c_xy = competition(trait,trait2,model)
        c_yy = competition(trait2,trait2,model)
        c_yx = competition(trait2,trait,model)
        if (c_xx*c_yy - c_xy*c_yx) == 0:
            raise ValueError('The traits x={} and y={} cannot coexsist on a long term.'.format(trait,trait2)) 
        else:
            return (r_x*c_yy - r_y*c_xy)/(c_xx*c_yy - c_xy*c_yx)
    
def invasion_fitness(y,x,model):
    """Calculates the invasion fitness of a trait y in an trait x monomorphic
    population at equilibrium size.
    
    """
    
    return birth_rate(y,model)-death_rate(y,model)-competition(y,x,model)*pop_equilibrium(x,model)
