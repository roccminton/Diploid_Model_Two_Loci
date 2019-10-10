"""Evolve a given population according to the standard model of adaptive 
dynamics, called BPDL (Boker, Pacala, Diekmann and Law), with or without mutation
and death due to natural reasons or due to competition.

The scripts expects that a model name is passed as an
argument. The model name must correspond to a file called
``[model_name].json`` and ``[model_name].py`` in the "IN_MODEL_SPECS" directory. 
In this way it is easy to extend the simulations to other model specifications if necessary.


"""

import json
import pickle
import importlib
import copy
import numpy as np
import os
from datetime import datetime
import time
import random

def setup_population():
    """Initializes the population by setting the time to zero and setup the *count* statistics which
    consits of a dictionary with six entries. Time, the trait statistic, which is itselfe
    a dictionary consisting of one entrie for every existing trait value and the corres-
    ponding numer of individuals within the current population carrying this trait. Then the same
    dictionary but with the total birthrates resp. deathrates as values. And last the
    overall population size and a mutation counter which increases by one whenever a mutation occurs. 

    """
    
    #Setup the initial population size which is handed in via the imported model. 
    #Every family starts with the same number of homogeneous members without diseses
    initial_population_size = int(model['K']*model['initial_rel_quantity_per_family'])*model['number_of_families']
    
    #Setup initial  trait_dic
    initial_trait_dic = {}
    initial_allele_config = ((0,)*model['number_of_loci'],)*2
    
    #loop over number of families and setup dictionary iteratively
    for family in range(1,model['number_of_families']+1):
        family_trait = (family,family,True)
        initial_trait_dic[family_trait] = [
                int(model['K']*model['initial_rel_quantity_per_family']),
                {initial_allele_config: int(model['K']*model['initial_rel_quantity_per_family'])}
                ]
    
    #setup the summerized statistics iteratively
    pop_state = {
        ('homo',0):initial_population_size,
            }

    #choose at random handed in amount of individuals to already carry the disease      
    for _ in range(int(initial_population_size*model['initial_rel_mutation_load'])):
        #choose a trait at random
        while True:
            choosen_family_trait = random.choice(list(initial_trait_dic.keys()))
            #if the randomly choosen family trait already is not able to reproduce choose one class of these non-reproductive
            #individuals and add one to the next class above
            if not choosen_family_trait[2]:
                choosen_category = random.choice(list(initial_trait_dic[choosen_family_trait][1].keys()))
                #reduce the number of individuals in this category by one
                _substract_and_check_if_empty(
                        dictionary = initial_trait_dic[choosen_family_trait][1],
                        key = choosen_category,
                        sub = 1
                        )
                #add the new individual to the dictionary in the category above the choosen one
                _check_if_key_exists_and_add(
                        dictionary = initial_trait_dic[choosen_family_trait][1],
                        key = choosen_category + 1,
                        add = 1
                        )
                #update the summerized statistics
                #add to new category
                _check_if_key_exists_and_add(
                    dictionary = pop_state,
                    key = ('homo',choosen_category + 1),
                    add = 1
                    )
                #substract from old category
                pop_state[('homo',choosen_category)] -= 1

                break
            
            choosen_allele_config = random.choice(list(initial_trait_dic[choosen_family_trait][1].keys()))
            choosen_chromosom = random.randint(0,1)
            choosen_locus = choose_locus(model['gene_distribution'])
            
            #check if the choosen allele at the chosen locus and chromosome didn't mutate yet
            if choosen_allele_config[choosen_chromosom][choosen_locus] == 0:
                #if so take away one of the individual from the initial list  
                #reduce the number of individuals in this category by one
                _substract_and_check_if_empty(
                        dictionary = initial_trait_dic[choosen_family_trait][1],
                        key = choosen_allele_config,
                        sub = 1
                        )
        
                #change the choosen allele configuration to carry the mutation at the randomly chosen site
                new_allele_config = _change_one_allele(
                        old_allele_config = choosen_allele_config,
                        chromosome = choosen_chromosom,
                        locus = choosen_locus,
                        new_type = 1)
                old_mutation_load = _ind_mut_load(choosen_allele_config)
                new_mutation_load = _ind_mut_load(new_allele_config)
                #update the summerized statistics 
                _check_if_key_exists_and_add(
                    dictionary = pop_state,
                    key = ('homo',new_mutation_load),
                    add = 1
                    )
                pop_state[('homo',old_mutation_load)] -= 1
                
                #check if the newly created type can still reproduce
                if new_allele_config[1-choosen_chromosom][choosen_locus] == 0:
                    #then add the ney type to the dictionary
                    _check_if_key_exists_and_add(
                            dictionary = initial_trait_dic[choosen_family_trait][1],
                            key = new_allele_config,
                            add = 1
                            )
                else:
                    #substract the individual also from the total family population
                    initial_trait_dic[choosen_family_trait][0] -= 1
                    #delete the family trait if extinct
                    if initial_trait_dic[choosen_family_trait][0] <= 0:
                        del initial_trait_dic[choosen_family_trait]
                    #create the new marked family trait
                    new_family_trait = choosen_family_trait[:2] + (False,)
                    #and add the new type marked as not reproductable to the dictionary
                    initial_trait_dic[new_family_trait] = [
                            1,
                            {new_mutation_load:1}
                            ]

                break
    
    #Setup the initial total birth and death rates per species with respective population sizes
    initial_birth_rates = {x_0: spc.birth_rate(x_0,model)*initial_trait_dic[x_0][0] for x_0 in initial_trait_dic.keys()}
    initial_death_rates = {x_0: spc.death_rate(x_0,model)*initial_trait_dic[x_0][0] for x_0 in initial_trait_dic.keys()}
    
    #Setup the initial compatibility per species including the additional death term for competition
    initial_compatibility_rates = {}

    for trait1 in initial_trait_dic.keys():
        initial_compatibility_rates[trait1] = {x_0: spc.birth_rate(x_0,model)*spc.rep_comp(trait1,x_0,model)*initial_trait_dic[x_0][0] for x_0 in initial_trait_dic.keys()}
        for trait2 in initial_trait_dic.keys():
            initial_death_rates[trait1] += (initial_trait_dic[trait1][0]/model['K'])*spc.competition(trait1,trait2,model)*initial_trait_dic[trait2][0]
                        
    #Initialize the initial statistic
    count = {
        'time':0,
        'trait': initial_trait_dic,
        'pop_state' : pop_state,
        'active_families': list(range(1,model['number_of_families']+1)),
        'mutation_counter':0,
        'mutation_load':int(initial_population_size*model['initial_rel_mutation_load']),
        'birth_rates': initial_birth_rates,
        'death_rates': initial_death_rates,
        'compatibility_rates': initial_compatibility_rates,
        }
    
    return count


def choose_species(rates):
    """Chooses a  trait at random according to its rates. The handed in variable 
    *rates* must be a dictionary with the traits to choose from as keys and the 
    probabilities to be chosen as values.
    
    """
    rnd = random.random() * sum(rates.values())
    for trait, rate in rates.items():
        rnd -= rate
        if rnd <= 0:
            return trait
        
def choose_locus(gene_distribution):
    """Chooses a gene section at random according to its sizes. The handed in variable
    *gene_distribution* must be a list of arbitrarily lenght which values sum up to
    the *number_of_loci* variable in the model dictionary.
    
    """
    if gene_distribution == 'unif':
        return random.choice(list(range(model['number_of_loci'])))
    elif type(gene_distribution) == list:
        rnd = np.random.random_sample() * sum(gene_distribution)
        for index, lenght in enumerate(gene_distribution):
            rnd -= lenght
            if rnd <= 0:
                return index
    else:
        raise ValueError('Unknown gene distribution: {}.'.format(gene_distribution))

    
def update_competition_and_compatibility(trait,event,manual_k=0):
    """Updates the competition pressure for every species when a single new individuum 
    of type *trait* enters/exits the population. Whether the individuum gets born or
    dies is specified in the event variable which can either be "birth" or "death"
    Note that the total competition pressure of a species x is given by 
        count(x)/K * sum(competition(x,y)*count(y) for y in X)
    Attention: The competition pressure of a species with itselfe in-/decreases quadratically,
    eg if count(x) -> count(x)+1 the competition pressure only induced by the species itselfe
    increases by competition(x,x)/K * 2(count(x)+1)-1
    """
    #Save global variables locally to speed up loop
    trait_list = dict([(x_0,count['trait'][x_0][0]) for x_0 in count['trait'].keys()])
    death_rates = count['death_rates']
    compatibility_rates = count['compatibility_rates']
    local_model = model
    
    #Set the scale parameter wether the competition pressure should be increased or reduced
    if event == 'birth':
        k = 1
    elif event == 'death':
        k = -1
    elif event == 'split':
        k = manual_k
    else:
        raise ValueError('The handed in event in *update_competition* should '/
                         +'be one of *birth* or *death*, but it is event='.format(event))
    
    #save factors for death rate updates:
    fac = k*(1/model['K'])
    
    for other_trait, size in trait_list.items():
        C = fac*size
        #updating competition preassure is not necessary for family splitting, since the new family will have
        #the same competition preasure from the other families and vice verca
        if event != 'split':
            #update death rate
            death_rates[trait] += C*competition(trait,other_trait,local_model)
            death_rates[other_trait] += C*competition(other_trait,trait,local_model)
        #update compatibility rate
        compatibility_rates[other_trait][trait] += k*spc.birth_rate(trait,local_model)*spc.rep_comp(other_trait,trait,local_model)
    
    if event != "split":    
        #since the species size for the handed in trait was updated before we need to substrainitial_birth_rates[x_0ct the overrun
        death_rates[trait] -= ((k**2)/model['K'])*competition(trait,trait,local_model)
    
    #apply the changes
    count['death_rates'] = death_rates
    count['compatibility_rates'] = compatibility_rates

def update_summerized_statistic(family_trait,mutation_load,k):
    """Updates the summerized statistics in which we only count the number of individuals with homogeneous/
    heterogeneous family labelling and of those these carrying the mutation resp. not.
    Here *k* must be or +1 or -1 indicating or the birth or the death of the individual with the handed in 
    trait. The trait must consist of a tuple with 4 entries. The first two giving the diploid family labelling,
    wheras the second two giving the diploid allele under mutation, where 0 means no mutation and 1 means the 
    individual carries the mutation.
    
    """    
    #if heterogeneous
    if family_trait[0] != family_trait[1]:
        #and not carrying the mutation
        _check_if_key_exists_and_add(
                dictionary = count['pop_state'],
                key = ('hetero',mutation_load),
                add = k
                )
    #if homogeneous
    elif family_trait[0] == family_trait[1]:
        #and not carrying the mutation
        _check_if_key_exists_and_add(
                dictionary = count['pop_state'],
                key = ('homo',mutation_load),
                add = k
                )
    else:
        if k == 1:
            event = 'birth'
        elif k == -1:
            event = 'death'
        else:
            raise ValueError('Handed in parameter k={} in update_summerized_statistic is invalid'.format(k))
        raise ValueError('The key {} in the update_summerized_statistic during a {} event is invalid.'.format(family_trait,event))

    
def competition(x,y,model):
    return spc.competition(x,y,model)

def execute_mutation(allele_config):
    """Executes a mutation at the handed in parent trait and returns the new trait.
    
    """
    #choose allele to mutate and execute mutation
    choose_chromosome = np.random.randint(2, size=1)[0]
    choosen_locus = choose_locus(model['gene_distribution'])
    allele_config = _change_one_allele(
            old_allele_config = allele_config,
            chromosome = choose_chromosome,
            locus = choosen_locus,
            new_type = spc.trait_mutation(allele_config[choose_chromosome][choosen_locus],model)
            )
    #update mutation counter
    count['mutation_counter'] += 1
    
    return allele_config

def intrachromosomal_recombination(allele_configuration):
    """Intrachromosomal recombination of the handed in allele_configuration means a random
    rearrangement of the chromosomes, which are specified via the model['chromosome_cuts'] and 
    model['gene_distribution'] lists.
    This function returns the haploid reshuffled allele_configuration. Notice that the handed in 
    allele_configuration is diploid!
    
    """
    #generate n_chromosome long sequence of 0 and 1 to choose from maternal resp. paternal chromosomes
    choice_sequence = list(np.random.randint(0,2,model['n_chromosome']))
    
    #setup new_allele_configuration
    new_allele_configuration = ()
    
    #go through each chromosome an add to the building up allele configuration or the
    #maternal or paternal chromosome, depending on the choice sequence
    for counter, chromosome_choice in enumerate(choice_sequence):
        new_allele_configuration += allele_configuration[chromosome_choice][model['chromosome_cuts'][counter]:model['chromosome_cuts'][counter+1]]
    
    return new_allele_configuration
    
    
def birth(birth_time, birth_rates):
    """ Returns the summerized statistics *count* after a new child was born at the handed in time
    where the individual trait value is chosen at random according to their birthrate.
    
    """
    #choose parent at random according to birth rates
    parent_family_trait = choose_species(birth_rates)
    
    #choose partner at random according to the compatibility rates
    while True:
        partner_family_trait = choose_species(count['compatibility_rates'][parent_family_trait])
        #check if partner is able to reproduce
        if partner_family_trait[2]:
            break     
    
    #random choice of parental family alleles
    choose_family_allele = np.random.randint(2, size=2)
    new_family_trait = tuple(sorted((parent_family_trait[choose_family_allele[0]],partner_family_trait[choose_family_allele[1]])))
    
    #choose parents and partners allele configuration according to their frequency 
    parent_allele_config = choose_species(count['trait'][parent_family_trait][1])
    partner_allele_config = choose_species(count['trait'][partner_family_trait][1])
    
    #if intrachromosomal recombination is turned on reshuffel the maternal and paternal chromosomes
    if model['intrachromosomal_rec'] == "True":
        #interchromosomal recombination between maternal and paternal information
        parent_allele_config = intrachromosomal_recombination(parent_allele_config)
        partner_allele_config = intrachromosomal_recombination(partner_allele_config)
        new_allele_config = tuple(sorted((parent_allele_config,partner_allele_config)))
    #if recombination is off, than just choose at random one of the two chromosome sets
    elif model['intrachromosomal_rec'] == "False":
        #random choice of parental allele combinations (without recombination)
        choose_allele_config = np.random.randint(2, size=2)
        new_allele_config = tuple(sorted((parent_allele_config[choose_allele_config[0]],partner_allele_config[choose_allele_config[1]])))
    else:
        raise ValueError('Unknown intrachromosomal recombination attribute: {}.'.format(model['intrachromosomal_rec']))
    
    #time is updated
    count['time'] += birth_time
    
    #if the mutation rate is measured per bais pair generate a poisson random variable and 
    #execute that many mutations, else just do one if the probability occures
    if model['inv_mutation_rate'] == 'per bp':
        N_mutation = np.random.poisson(model['mutation_rate'])
        for _ in range(N_mutation):
            new_allele_config = execute_mutation(new_allele_config)
    else:
        #calculate mutation probabilities
        p = spc.trait_mutation_probability(model)
        U = np.random.random_sample()
        #trait mutation, hence U in (0, p)
        if U <= p:
            new_allele_config = execute_mutation(new_allele_config)
        
    #update mutation load
    mutation_load = _ind_mut_load(new_allele_config)
    count['mutation_load'] += mutation_load
    
    #check if newly created allele configuration is able to reproduce an set the boolian accordingly
    reproduction = True
    for a,A in zip(new_allele_config[0],new_allele_config[1]):
        if a >= 1 and A >= 1:
            reproduction = False
            new_allele_config = mutation_load
            if mutation_load < 2:
                raise ValueError('Individuum {} was marked not reproducible with a mutation load of {}.'.format(new_family_trait,mutation_load))
            break 
    #save the reproduction trait in the family trait
    new_family_trait += (reproduction,)
    
    #add new trait to trait dict and update birth, death, and compatibility rates
    if new_family_trait in count['trait'].keys():
        count['trait'][new_family_trait][0] += 1
        _check_if_key_exists_and_add(
                dictionary = count['trait'][new_family_trait][1],
                key = new_allele_config,
                add = 1
                )
        count['birth_rates'][new_family_trait] += spc.birth_rate(new_family_trait,model)
        count['death_rates'][new_family_trait] += spc.death_rate(new_family_trait,model)    
    else:
        count['trait'][new_family_trait] = [
                1,
                {new_allele_config:1}
                ]
        count['birth_rates'][new_family_trait] = spc.birth_rate(new_family_trait,model)
        count['death_rates'][new_family_trait] = spc.death_rate(new_family_trait,model)
        #setup a new dictionary entry for the newly arrived trait
        count['compatibility_rates'][new_family_trait] = {}
        for other_trait in count['trait'].keys():
            #setup the new register in the compatibility matrix, but substract one factor, since it will be updated again in the next step
            count['compatibility_rates'][new_family_trait][other_trait] = spc.birth_rate(other_trait,model)*spc.rep_comp(new_family_trait,other_trait,model)*(count['trait'][other_trait][0]-1)
            count['compatibility_rates'][other_trait][new_family_trait] = 0
        
    #update summerized statistics
    update_summerized_statistic(new_family_trait,mutation_load,1)    
    
    #Update competition pressure and reproductive compatibility in whole population
    #needs to be done after increasing the species size
    update_competition_and_compatibility(new_family_trait,'birth')

def death(death_time, death_rates):
    """ Returns the summerized statistics *count* after a individuum died at the handed in time
    where the individual trait value is chosen at random according to their deathrate.
    
    """
    #choose fey
    fey_family_trait = choose_species(death_rates)   
    #choose fey allele configuration
    fey_allele_config = choose_species(count['trait'][fey_family_trait][1])
    
    #mutation load update
    if fey_family_trait[2]:
        mutation_load = _ind_mut_load(fey_allele_config)
    else:
        mutation_load = fey_allele_config
    count['mutation_load'] -= mutation_load
    
    #total statistic is updated
    count['trait'][fey_family_trait][0] -= 1
    _substract_and_check_if_empty(
            dictionary = count['trait'][fey_family_trait][1],
            key = fey_allele_config,
            sub = 1
            )

    #summerized statistics is updated
    update_summerized_statistic(fey_family_trait,mutation_load,-1)
    
    #time is updated
    count['time'] += death_time

    #birth and death rates are reduced
    count['birth_rates'][fey_family_trait] -= spc.birth_rate(fey_family_trait,model)
    count['death_rates'][fey_family_trait] -= spc.death_rate(fey_family_trait,model)
    
    #Update competition pressure in whole population
    #needs to be done after decreasing the species size
    update_competition_and_compatibility(fey_family_trait,'death')
        
    #erase extinct species from dictionary
    if count['trait'][fey_family_trait][0] == 0:
        del count['trait'][fey_family_trait]
        del count['birth_rates'][fey_family_trait]
        del count['death_rates'][fey_family_trait]
        del count['compatibility_rates'][fey_family_trait]
        for other_trait in count['trait'].keys():
            del count['compatibility_rates'][other_trait][fey_family_trait]
            
    #check if there are still members of the family left
    family_a = fey_family_trait[0]
    family_b = fey_family_trait[1]
    found_a = False
    found_b = False
    
    #loop over all active traits and check if the family lable is still represented
    for lable in count['trait']:
        if family_a in lable[:2]:
            found_a = True
        if family_b in lable[:2]:
            found_b = True
        if found_a == found_b == True:
            break
        
    #if the checking loop got through without breaking remove the extinct 
    #family from the list of active families
    if found_a == False:
        count['active_families'].remove(family_a)
    if found_b == False and (family_a != family_b):
        count['active_families'].remove(family_b)

def next_event_time():
    """ Returns a sample of an exponential distributed random variable representing the minimum of all
    the independent exponential distributed birth- and death times. Remark that the minimum over indep.
    exponential distributed rv. is again an exponential distributed rv. with parameter the sum of the 
    initial parameters of the exp. rvs.

    """

    birth_rates = count['birth_rates']
    death_rates = count['death_rates']
    
    #Draw an exponential distribiuted random variable with the sum of birth resp. deaht rates as parameter
    #if the all birth or death rates happen to be zero setz the first birth resp. death time to infinity
    #if both, all birth AND all death rates are zero somehow the evolution stopped and return a break command to the main loop
    try:
        first_birth_time = np.random.exponential(sum(birth_rates.values())**(-1))
    except ZeroDivisionError:
        first_birth_time = np.inf
    try:
        first_death_time = np.random.exponential(sum(death_rates.values())**(-1))
    except ZeroDivisionError:
        first_death_time = np.inf
        
    if first_birth_time == np.inf and first_death_time == np.inf:
        return first_birth_time, 0, 'stop'
    elif first_birth_time < first_death_time:
        return first_birth_time, birth_rates, 'birth'
    else:
        return first_death_time, death_rates, 'death'

def one_step_evolution():
    """ Executes one demographical event and returns the new population.

    """

    time, rates, event = next_event_time()

    if event == 'birth':
        birth(time, rates)
    elif event == 'death':
        death(time, rates)
    elif event == 'stop':
        count['stop'] = True
    else: 
        raise ValueError(
            'There is something wrong in the "next_event" function in src.analysis.evolve.'
            )

def display_time(seconds, granularity=2):
    """Converts n seconds into w weeks, d days, h hours, m minutes and s seconds.
    
    """
    
    result = []
    
    intervals = (
    ('weeks', 604800),  # 60 * 60 * 24 * 7
    ('days', 86400),    # 60 * 60 * 24
    ('hours', 3600),    # 60 * 60
    ('minutes', 60),
    ('seconds', 1),
    )

    for name, number in intervals:
        value = seconds // number
        if value:
            seconds -= value * number
            if value == 1:
                name = name.rstrip('s')
            result.append("{} {}".format(value, name))
    return ', '.join(result[:granularity])

def extract_summerized_statistics():
    """Reduce the count dictionary to the keys and values one wants to keep for 
    saving and plotting.
    
    """
    fam_distribution = {}
    #build up abs frequency list
    for trait in count['trait'].keys():
        _check_if_key_exists_and_add(
                dictionary = fam_distribution,
                key = count['trait'][trait][0],
                add = 1
                )
    
    summerized_statistics = {
            'time': count['time'],
            'pop_state' : count['pop_state'],
            'active_families': len(count['active_families']),
            'mutation_counter': count['mutation_counter'],
            'mutation_load': count['mutation_load'],
            'number_false' : sum([count['trait'][family_trait][0] for family_trait in filter(lambda key: key[2]==False, count['trait'].keys())]),
            'number_of_types':len(count['trait']),
            'fam_distribution': fam_distribution
            }
    
    return summerized_statistics

def build_up_averaged_statistics(averaged_statistics,summerized_statistics,num_runs):
    """Let the simulation run for *num_runs* times this function averages the extracted
    summerized statistics into an averaged version of it.
    
    """

    for key in averaged_statistics.keys():
        #if the values are just numbers add the averaged stats to the corresponding key
        if key not in ['fam_distribution','pop_state']:
            averaged_statistics[key] += summerized_statistics[key]/num_runs
        #else if the values are dictionaries again loop over every key and add the single averaged stats
        else:
            for keyy in summerized_statistics[key].keys():
                #check if the key in fam_distribution or pop_state already exists in the averaged version
                if keyy in averaged_statistics[key].keys():
                    #if yes just add the averaged quantity
                    averaged_statistics[key][keyy] += summerized_statistics[key][keyy]/num_runs
                else:
                    #if not create a new entry with the averaged stats
                    averaged_statistics[key][keyy] = summerized_statistics[key][keyy]/num_runs
        
        
    return averaged_statistics

def _check_if_key_exists_and_add(dictionary,key,add):
    """Checks if the handed in key exists in the handed in dictionary and adds
    the handed in amount to the key if possible or creates a new key in the dictonary
    with the handed in amount as value
    
    """
    if key in dictionary.keys():
        dictionary[key] += add
    else:
        dictionary[key] = add
        
def _substract_and_check_if_empty(dictionary,key,sub):
    """The handed in key must be a indeed a key of the dictionary. We substract the handed
    in amount *sub* from the value of the dictionaries key and check afterwards if the 
    value is still positive. If not we delete the entry from the dictionary.
    
    """
    
    dictionary[key] -= sub
    
    if dictionary[key] <= 0:
        del dictionary[key]

def _change_one_allele(old_allele_config,chromosome,locus,new_type):
    """Changes the Allele Configuration at the handed in chromosome and locus to *new_allele* and
    returns the newly created configuration. The allele configuration should be a tuple consisting 
    of two entries. These two entries should be tuples as well of the same size.

    """
    new_chromosome = old_allele_config[chromosome][:locus] + (new_type,) + old_allele_config[chromosome][locus+1:]
    
    return tuple(sorted((new_chromosome,old_allele_config[1-chromosome])))
    
def _ind_mut_load(allele_config):
    """Calculates the individuals mutation load from the handed in allele configuration.
    
    
    """
    #convert the tuple of tuples -> tuple -> sum
    return sum(sum(allele_config, ()))

def K_change(time,k,k_changes):
    """Gives the carrying capacity with respect to the time. 
    
    """
    #check if list is empty 
    if not model['K_change']:
        return k, k_changes
    elif time >= model['K_change'][-1]:
        return k, k_changes
    elif time + model['step'] >= model['K_change'][k_changes] and time < model['K_change'][k_changes]:
        k_changes += 1
        k = model['K_list'][k_changes]
        model['number_of_families'] = model['number_of_families_list'][k_changes]
        return k, k_changes
    else:
        return k, k_changes
    
def _chromosome_cuts(gene_distribution,n_chromosome):
    """Returns a list of where to cut the gene distribution list to get *n_chromosomes* of
    almost even length. Go through the gene distribution list and every time the accumulated
    sums of the entries exceed the average chromosome length the list index is saved in the 
    newly created *cuts* list. 
    We add 0 in the beginning and len(gene_distribution)+1 in the end of the *cuts* list so
    that chromosomes can get recovered via slicing.
    Note that all but the last chromosome will be bigger than the average chromosome size.
    
    """
    average_chromosome_length = sum(gene_distribution)/n_chromosome

    run = 0
    cuts = [0]
    
    for counter, value in enumerate(gene_distribution):
        run += value
        if run >= average_chromosome_length:
            cuts.append(counter+1)
            run = 0
    
    cuts.append(len(gene_distribution)+1)
    
    if len(cuts) != n_chromosome + 1:
        raise KeyError("Chromosome cuts went wrong. Got {} chromosomes instead of {}.".format(len(cuts)-1,n_chromosome))
    
    return cuts

def _divide_dictionary_into_two_at_random(dictionary,total_sum=None):
    """Iterates over all the values in the dictionary and returns the index of the entry
    where the sum of the values before and after the index are almost equal. Additionaly 
    return the sum of the first half.
    
    """
    if total_sum == None:
        total_sum = sum(dictionary.values())
    elif total_sum == sum(dictionary.values()):
        #set the approx size of each halfe of the dictionaries
        total_sum = total_sum // 2
        #safe the list of keys of the dictionary
        keys = list(dictionary.keys())
        #and reshuffle them in a random order
        random.shuffle(keys)
        #make a copy of the dictionary and modify both simultaniously
        second_dict = copy.deepcopy(dictionary)
        #iterate over every key and split the dictionaries a part
        for key in keys:
            if total_sum > 0:
                size = np.random.binomial(min(total_sum,dictionary[key]),1/2)
#                if dictionary[key] == 1:
#                    size = np.random.binomial(min(total_sum,dictionary[key]),1/2)
#                else: 
#                    size = dictionary[key] // 2
                second_size = dictionary[key]-size
                _substract_and_check_if_empty(
                        dictionary=dictionary,
                        key=key,
                        sub=size)
                _substract_and_check_if_empty(
                        dictionary=second_dict,
                        key=key,
                        sub=second_size)
                total_sum -= size
            else:
                del dictionary[key]
        return dictionary, second_dict     
    else:
        raise ValueError('Handed in sum of the dictionaries values ({}) does not match the real sum ({})'.format(
                total_sum,
                sum(dictionary.values())))
    

def check_family_sizes():
    """Runs over every family entry, checks the number of individuals and if the families
    become to large split them up in new, homogene families.
    
    """

    #set up the maximum family size
    max_fam_size = 2*model['K']/model['number_of_families']

    #save the count dictionary to evoid Runtime Error because of changeing dictionary sizes during iteration
    for trait in list(count['trait'].keys()):
        if count['trait'][trait][0] >= max_fam_size:
            #delete old family from summerized statistics
            mutation_load_dict = _extract_mutation_load_dict(trait_dict=count['trait'][trait][1])
            old_birth_rate = count['birth_rates'][trait]
            old_death_rate = count['death_rates'][trait]
            for mutation_load, size in mutation_load_dict.items():
                update_summerized_statistic(trait,mutation_load,-size)
            
            #safe old family size
            old_family_size = count['trait'][trait][0]
            
            #create two new dictionaries for the new types with half the enties each 
            old_key_dict, new_key_dict = _divide_dictionary_into_two_at_random(count['trait'][trait][1],total_sum=old_family_size)
            
            #set the new family trait
            new_trait = (model['latest_family_index']+1,model['latest_family_index']+1,trait[2])
            
            #safe sizes of new families
            old_part = sum(old_key_dict.values())
            new_part = sum(new_key_dict.values())
            
            #check if one family size is zero
            if old_part == 0 or new_part == 0:
                raise ValueError('Invalid splitting of family {} into {} and {}.'.format(trait,old_part,new_part))
           
            #check if total size hasn't changed
            if old_part + new_part != old_family_size:
                raise ValueError('Somebody was left behind or countet twice during split of family. Population size before {} after {} split.'.format(old_family_size,old_part + new_part))
            
            #save both new types 
            count['trait'][trait] = [old_part, old_key_dict]
            count['trait'][new_trait] = [new_part, new_key_dict]
            
            total_birth_rate = sum(count['birth_rates'].values())
            total_death_rate = sum(count['death_rates'].values())
            
            #update birth and deathrates 
            #notice that competition rates for the other families does not have to be updated
            #since the old and the new family execute the same competition preassure on the other individuals
            count['birth_rates'][trait] = old_birth_rate*(old_part/old_family_size)
            count['birth_rates'][new_trait] = old_birth_rate*(new_part/old_family_size)
            count['death_rates'][trait] = old_death_rate*(old_part/old_family_size)
            count['death_rates'][new_trait] = old_death_rate*(new_part/old_family_size)
            
            #setup a new compatibility rates dictionary for new_trait with the old rates
            count['compatibility_rates'][new_trait] = {x_0: spc.birth_rate(x_0,model)*spc.rep_comp(new_trait,x_0,model)*count['trait'][x_0][0] for x_0 in count['trait'].keys()}
            for other_trait_dict in count['compatibility_rates'].values():
                other_trait_dict[new_trait] = 0
            
            #update compatibility rates (update of the competition is switched off for the split event)
            update_competition_and_compatibility(trait,"split",manual_k=(-1)*new_part)
            update_competition_and_compatibility(new_trait,"split",manual_k=new_part)
            
            #check if birth and death rates have been updated correctly
            if round(total_birth_rate,2) != round(sum(count['birth_rates'].values()),2):
                raise ValueError('Total birthrate has changed from {} to {}.'.format(total_birth_rate,sum(count['birth_rates'].values())))
            if round(total_death_rate,2) != round(sum(count['death_rates'].values()),2):
                raise ValueError('Total death rate has changed from {} to {}.'.format(total_death_rate,sum(count['death_rates'].values())))

            #extract mutation load dictionaries
            old_mutation_load_dict = _extract_mutation_load_dict(trait_dict=old_key_dict)
            new_mutation_load_dict = _extract_mutation_load_dict(trait_dict=new_key_dict)
            
            #update summerized statistics
            for mutation_load, size in old_mutation_load_dict.items():
                update_summerized_statistic(trait,mutation_load,size)
            for mutation_load, size in new_mutation_load_dict.items():
                update_summerized_statistic(new_trait,mutation_load,size)

            #update latest family index
            model['latest_family_index'] += 1
            #update list of active families
            count['active_families'].append(model['latest_family_index'])
            
            
            

def _extract_mutation_load_dict(trait=None,trait_dict=None):
    """Extracts the mutatio load dictionary from the dictionary *count['trait']* with the handed
    in key *trait* or from a dictionary *trait_dictionary* of the same form. The dictionary must
    have as keys the allele_configuration from which we can extract the mutation load and as value
    the number of individuals carrying this exact configuration. 
    The output is a dictionary with keys the mutation load and as values the frequency of the mutation
    load in the handed in key.
    
    """
    if trait != None:
        dictionary = count['trait'][trait][1]
    elif trait_dict != None:
        dictionary = trait_dict
    else:
        raise ValueError("_extract_mutation_load_dict must be raised with at least one argument, None was given.")
    
    mutation_load_dict = {}
    
    #iterate over every allele configuration
    for configuration, size in dictionary.items():
        mutation_load = _ind_mut_load(configuration)
        _check_if_key_exists_and_add(dictionary=mutation_load_dict,
                                     key=mutation_load,
                                     add=size
                                     )
    #check if sizes fit
    if sum(dictionary.values()) != sum(mutation_load_dict.values()):
        raise ValueError("There ist a mistake in the _extract_mutation_load_function. {] nicht gleich {}.".format(sum(dictionary.values()),sum(mutation_load_dict.values())))
    
    return mutation_load_dict

def debug_population_size():
    """Checks if the population size is updated correctly.
    
    """
    sum_over_traits = sum(count['trait'][x_0][0] for x_0 in count['trait'].keys())
    sum_over_pop_state = sum(count['pop_state'].values())
    
    if sum_over_traits != sum_over_pop_state:
        raise ValueError('There is a difference in the population sizes. Sum over all traits gives {}, whereas sum over the population state gives {}'.format(sum_over_traits,sum_over_pop_state))

def debug_trait_dict():
    """Iterates over every trait dict an checks if the sizes fit.
    
    """
    for trait, trait_dict in count['trait'].items():
        dict_sum = sum(trait_dict[1].values())
        if dict_sum != trait_dict[0]:
            raise ValueError('The trait {} was updated incorrectly. Shown size is {}, but actual size is {}.'.format(trait,trait_dict[0],dict_sum) )

if __name__ == "__main__":
    model_name = "Families"

    
    #Absolute dir the script is in. Needs to be changed manually if the file is moveed
    script_dir = r'/home/larocca/Dokumente/Simulationen/Diploid_Model_Two_Loci'
    
    #load model specific variables such as carrying capacity, initial trait and marker, etc.
    rel_path = "model_specifications/{}.json".format(model_name)
    abs_file_path = os.path.join(script_dir, rel_path)
    model = json.load(open(abs_file_path), encoding="utf-8")
    
    #import model specific functions such as birthrate, mutation operator and probabilities, etc.
    spc = importlib.import_module('model_specifications.{}'.format(model_name))
    
    #set clock to zero
    t_start = time.clock()
    
    #set global variables
    model['K'], K_changes = model['K_list'][0], 0
    model['number_of_families'] = model['number_of_families_list'][0]
    model['latest_family_index'] = model['number_of_families']
    
    #check if runtime is big enough for changes in K
    if len(model['K_change']) > 0:
        if model['K_change'][-1] > model['t_max']:
            raise ValueError("The last changing time for K: {} is bigger than the total runtime: {}.".format([model['K_change'][-1]],model['t_max']))
    
    #setup the gene distribution if necessary
    if type(model['gene_distribution']) == list:
        model['gene_distribution'] = list(np.random.randint(model['gene_distribution'][0],model['gene_distribution'][1]+1,model['gene_distribution'][2]))
        model['number_of_loci'] = len(model['gene_distribution'])
        model['chromosome_cuts'] = _chromosome_cuts(model['gene_distribution'],model['n_chromosome'])
    #setup the mutation rate from the number of loci if necessary
    if model['inv_mutation_rate'] == "per bp":
        model['mutation_rate'] = sum(model['gene_distribution'])*model['mutation_rate']
    #setup recombination rate from the number of loci if necessary
    if "rec_rate" in model.keys():
        model['rec_rate'] = sum(model['gene_distribution'])*model['rec_rate']
        
    #setup averaged stats per time list with empty averaged stats in it
    empty_stats = {
        'time': 0,
        'pop_state' : {},
        'active_families': 0,
        'mutation_counter': 0,
        'mutation_load': 0,
        'number_false' : 0,
        'number_of_types':0,
        'fam_distribution': {}
        }
    
    #build up empty averaged_stats_per_time liste iteratively in for loop
    #to have an independent deepcopy of empty_stats dictionary as list entries
    averaged_stats_per_time = []
    for _ in range(model['t_max']+1):
        averaged_stats_per_time.append(copy.deepcopy(empty_stats))
        
    #set the variable for a family check
    family_check = 50
    
    run = 1
    #Loop over the number of runs one wants to average over
    while run <= model['num_runs']:        
        #check if the run terminated succsessfully and ran to the end without the population gettin extinct before
        terminated_successfully = True

        #set run clock to zero
        t_start_run = time.clock()
        
        # Load initial population and setup random times
        count = setup_population()
        
        #compromized statistic that stores time and number of traits and makers
        stats_per_time = [copy.deepcopy(extract_summerized_statistics())]
                
        # Run the main analysis up to some fixed time t_max and safe the statistics every *step* times
        for t in range(0,model['t_max'],model['step']):
            model['K'], K_changes = K_change(t,model['K'],K_changes)
            while count['time'] <= t:
                one_step_evolution()
                if sum(count['pop_state'].values())==0:
                    print('Population got extinct at time {}.'.format(count['time']))
                    terminated_successfully = False
                    break
                elif 'stop' in count.keys():
                    print('Evolution stopped at time {} with a total of {} individuals.'.format(count['time'],sum(count['pop_state'].values())))
                    terminated_successfully = False
                    break
            if sum(count['pop_state'].values())==0:
                terminated_successfully = False
                break
            if min(count['pop_state'].values()) < 0:
                print('One population state got negative')
                terminated_successfully = False
                break
            #every *family_check* steps check the family sizes and split them up eventually
            #only if families matter to speed up runtime
            if model['mating_scheme'][0] != "random":
                if t % family_check == 0:
                    check_family_sizes()
                    debug_trait_dict()
                    debug_population_size()
            #extract summerized statistics to safe in file
            stats_per_time.append(copy.deepcopy(extract_summerized_statistics()))
            print('Process: {}% of run {}/{} at {} (X_t = {} for t={})'.format(t*100/model['t_max'],run,model['num_runs'],datetime.time(datetime.now()),sum(count['pop_state'].values()),t))

        #get time difference for single run 
        t_end_run = time.clock()
        dt_run = t_end_run - t_start_run
          
        print('Evolution terminated successfully. At the final time t={} there were {} individuals alive. The simulation ran for {}.'.format(count['time'],sum(count['pop_state'].values()),display_time(dt_run)))
    
        if terminated_successfully:
            #average over the runs and safe it in new list
            for counter, (averaged_stats, summerized_stats) in enumerate(zip(averaged_stats_per_time,stats_per_time)):
                averaged_stats_per_time[counter] = copy.deepcopy(build_up_averaged_statistics(averaged_stats,summerized_stats,model['num_runs']))
            run += 1
        
        #don't reset the global variables in the last run
        if run < model['num_runs']:
            #reset global variables
            model['K'], K_changes = model['K_list'][0], 0
    
            
    #get time difference for total simulation 
    t_end = time.clock()
    dt = t_end - t_start
    
    print('The total {}/{} simulations ran for {}.'.format(model['num_runs'],model['num_runs'],display_time(dt)))
    print('The following parameters have been used:')
    print(model)
    
    # Store list with locations after each round in a pickle file
    rel_path = "out_analysis/data/statistic_{}(K={},nloci={},mating={},n_runs={}).pickle".format(
            model_name,
            model['K_list'],
            model['number_of_loci'],
            model['mating_scheme'],
            model['num_runs']
            )
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "wb") as out_file:
        pickle.dump(averaged_stats_per_time, out_file)
