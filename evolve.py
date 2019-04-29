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
    initial_trait = {}
    initial_state = []
    initial_population_size = int(model['K']*model['initial_rel_quantity_per_family'])*model['number_of_families']

    for family in range(1,model['number_of_families']+1):
        trait = tuple([family,family,0,0])
        initial_state.append(trait)
        initial_trait[trait] = int(model['K']*model['initial_rel_quantity_per_family'])
    
    #choose at random handed in amount of individuals to already carry the disease      
    for _ in range(int(initial_population_size*model['initial_rel_quantity_of_disease'])):
        #choose a trait at random
        while True:
            choosen_trait = random.choice(list(initial_trait.keys()))
            #check if the choosen trait is still represented in the population
            #and if it is a (0,0) trait
            if initial_trait[choosen_trait] > 0 and (choosen_trait[2]+choosen_trait[3]) == 0:
                break
        #if so take away one of the individual from the (0,0) list   
        initial_trait[choosen_trait] -= 1
        #and erase extinct species eventually from list
        if initial_trait[choosen_trait] == 0:
            del initial_trait[choosen_trait]
            initial_state.remove(choosen_trait)

        #and add the new (0,1) individual to it
        new_trait = tuple([choosen_trait[0],choosen_trait[1],0,1])
        if new_trait not in initial_state:
            initial_trait[new_trait] = 1
            initial_state.append(new_trait)
        else:
            initial_trait[new_trait] +=1
    
    #Setup the initial total birth rates per species with respective population sizes
    initial_birth_rates = {x_0: spc.birth_rate(x_0,model)*initial_trait[x_0] for x_0 in initial_state}
    
    #Setup the initial totla death rates and initial compatibility per species including competition
    initial_death_rates = {}
    initial_compatibility_rates = {}

    for trait1 in initial_state:
        initial_death_rates[trait1] = spc.death_rate(trait1,model)*initial_trait[trait1]
        initial_compatibility_rates[trait1] = {x_0: spc.birth_rate(x_0,model)*spc.rep_comp(trait1,x_0,model)*initial_trait[x_0] for x_0 in initial_state}
        for trait2 in initial_state:
            initial_death_rates[trait1] += (initial_trait[trait1]/model['K'])*spc.competition(trait1,trait2,model)*initial_trait[trait2]
            
    pop_state = {
            ('homo',0):int(initial_population_size*(1-model['initial_rel_quantity_of_disease'])),
            ('homo',1):int(initial_population_size*model['initial_rel_quantity_of_disease']),
            ('hetero',0):0,
            ('hetero',1):0,
                }
            
    #Initialize the initial statistic
    count = {
        'time':0,
        'trait': initial_trait,
        'pop_state' : pop_state,
        'active_families': list(range(1,model['number_of_families']+1)),
        'mutation_counter':0,
        'birth_rates': initial_birth_rates,
        'death_rates': initial_death_rates,
        'compatibility_rates': initial_compatibility_rates,
        'trait_list': initial_state
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
    
def update_competition_and_compatibility(trait,count,event):
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
    trait_list = count['trait']
    death_rates = count['death_rates']
    compatibility_rates = count['compatibility_rates']
    local_model = model
    
    
    #Set the scale parameter wether the competition pressure should be increased or reduced
    if event == 'birth':
        k = 1
    elif event == 'death':
        k = -1
    else:
        raise ValueError('The handed in event in *update_competition* should '/
                         +'be one of *birth* or *death*, but it is event='.format(event))
    
    #save factors for death rate updates:
    fac = k*(1/local_model['K'])
    
    for other_trait, size in trait_list.items():
        C = fac*size
        #update death rate
        death_rates[trait] += C*competition(trait,other_trait,local_model)
        death_rates[other_trait] += C*competition(other_trait,trait,local_model)
        #update compatibility rate
        compatibility_rates[other_trait][trait] += k*spc.birth_rate(trait,local_model)*spc.rep_comp(other_trait,trait,local_model)

    #since the species size for the handed in trait was updated before we need to substrainitial_birth_rates[x_0ct the overrun
    death_rates[trait] -= competition(trait,trait,local_model)/local_model['K']
    
    #apply the changes
    count['death_rates'] = death_rates
    count['compatibility_rates'] = compatibility_rates

def update_summerized_statistic(trait,count,k):
    """Updates the summerized statistics in which we only count the number of individuals with homogeneous/
    heterogeneous family labelling and of those these carrying the mutation resp. not.
    Here *k* must be or +1 or -1 indicating or the birth or the death of the individual with the handed in 
    trait. The trait must consist of a tuple with 4 entries. The first two giving the diploid family labelling,
    wheras the second two giving the diploid allele under mutation, where 0 means no mutation and 1 means the 
    individual carries the mutation.
    
    """
    #if heterogeneous
    if trait[0] != trait[1]:
        #and not carrying the mutation
        if (trait[2]+trait[3]) == 0:
            count['pop_state'][('hetero',0)] += k
        else:
            count['pop_state'][('hetero',1)] += k

    #if homogeneous
    elif trait[0] == trait[1]:
        #and not carrying the mutation
        if (trait[2] + trait[3]) == 0:
            count['pop_state'][('homo',0)] += k
        else:
            count['pop_state'][('homo',1)] += k
    else:
        if k == 1:
            event = 'birth'
        elif k == -1:
            event = 'death'
        else:
            raise ValueError('Handed in parameter k={} in update_summerized_statistic is invalid'.format(k))
        raise ValueError('The key {} in the update_summerized_statistic during a {} event is invalid.'.format(trait,event))

    
def competition(x,y,model):
    return spc.competition(x,y,model)
    
def birth(birth_time, birth_rates, count):
    """ Returns the summerized statistics *count* after a new child was born at the handed in time
    where the individual trait value is chosen at random according to their birthrate.
    
    """
    #choose parent at random according to birth rates
    parent_trait = choose_species(birth_rates)
    
    #choose partner at random according to the compatibility rates
    partner_trait = choose_species(count['compatibility_rates'][parent_trait])
    
    #random recombination of parental alleles
    choose_family_allele = np.random.randint(2, size=2)
    choose_disease_allele = np.random.randint(low=2, high=4, size=2)
    new_family_trait = tuple(sorted([parent_trait[choose_family_allele[0]],partner_trait[choose_family_allele[1]]]))
    new_disease_trait = tuple(sorted([parent_trait[choose_disease_allele[0]],partner_trait[choose_disease_allele[1]]]))
    
    #time is updated
    count['time'] += birth_time

    #calculate mutation probabilities
    p = spc.trait_mutation_probability(model)
    U = np.random.random_sample()

    #trait mutation, hence U in (0, p)
    if U <= p:
        #choose allele to mutate
        choose_allele = np.random.randint(2, size=1)[0]
        allele_to_mutate = new_disease_trait[choose_allele]
        #execute mutation on chosen allele
        new_allele = spc.trait_mutation(allele_to_mutate,model)
        #combine the parental allele with the mutated one
        new_disease_trait =  tuple(sorted([new_allele,new_disease_trait[1-choose_allele]]))

        #update mutation counter
        count['mutation_counter'] += 1
        
    new_trait = new_family_trait + new_disease_trait
        
    #add new trait to trait_list and update birth, death, and compatibility rates
    if new_trait not in count['trait_list']:
        count['trait_list'].append(new_trait)
    if new_trait in count['trait'].keys():
        count['trait'][new_trait] += 1 
        count['birth_rates'][new_trait] += spc.birth_rate(new_trait,model)
        count['death_rates'][new_trait] += spc.death_rate(new_trait,model)    
    else:
        count['trait'][new_trait] = 1
        count['birth_rates'][new_trait] = spc.birth_rate(new_trait,model)
        count['death_rates'][new_trait] = spc.death_rate(new_trait,model)
        #setup a new dictionary entry for the newly arrived trait
        count['compatibility_rates'][new_trait] = {}
        for other_trait in count['trait'].keys():
            #setup the new register in the compatibility matrix, but substract one factor, since it will be updated again in the next step
            count['compatibility_rates'][new_trait][other_trait] = spc.birth_rate(other_trait,model)*spc.rep_comp(new_trait,other_trait,model)*(count['trait'][other_trait]-1)
            count['compatibility_rates'][other_trait][new_trait] = 0
    
    #update summerized statistics
    update_summerized_statistic(new_trait,count,1)    
    
    #Update competition pressure and reproductive compatibility in whole population
    #needs to be done after increasing the species size
    update_competition_and_compatibility(new_trait,count,'birth')

def death(death_time, death_rates, count):
    """ Returns the summerized statistics *count* after a individuum died at the handed in time
    where the individual trait value is chosen at random according to their deathrate.
    
    """
    #choose fey
    fey_trait = choose_species(death_rates)
    
    #summerized statistic is updated
    count['trait'][fey_trait] -= 1
    
    #time is updated
    count['time'] += death_time

    #update summerized statistics
    update_summerized_statistic(fey_trait,count,-1)
    
    #birth and death rates are reduced
    count['birth_rates'][fey_trait] -= spc.birth_rate(fey_trait,model)
    count['death_rates'][fey_trait] -= spc.death_rate(fey_trait,model)
    
    #Update competition pressure in whole population
    #needs to be done after decreasing the species size
    update_competition_and_compatibility(fey_trait,count,'death')
        
    #erase extinct species from dictionary
    if count['trait'][fey_trait] == 0:
        del count['trait'][fey_trait]
        del count['birth_rates'][fey_trait]
        del count['death_rates'][fey_trait]
        del count['compatibility_rates'][fey_trait]
        for other_trait in count['trait'].keys():
            del count['compatibility_rates'][other_trait][fey_trait]
            
    #check if there are still members of the family left
    family_a = fey_trait[0]
    family_b = fey_trait[1]
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

def next_event_time(count):
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

def one_step_evolution(count):
    """ Executes one demographical event and returns the new population.

    """

    time, rates, event = next_event_time(count)

    if event == 'birth':
        birth(time, rates, count)
    elif event == 'death':
        death(time, rates, count)
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

    for name, count in intervals:
        value = seconds // count
        if value:
            seconds -= value * count
            if value == 1:
                name = name.rstrip('s')
            result.append("{} {}".format(value, name))
    return ', '.join(result[:granularity])

def extract_summerized_statistics(count):
    """Reduce the count dictionary to the keys and values one wants to keep for 
    saving and plotting.
    
    """
    summerized_statistics = {
            'time': count['time'],
            'pop_state' : count['pop_state'],
            'active_families': count['active_families'],
            'mutation_counter': count['mutation_counter'],
            }
    
    return summerized_statistics


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
            
    # Load initial population and setup random times
    count = setup_population()
    
    #compromized statistic that stores time and number of traits and makers
    stats_per_time = [copy.deepcopy(extract_summerized_statistics(count))]

    # Run the main analysis up to some fixed time t_max and safe the statistics every *step* times
    for t in range(0,model['t_max'],model['step']):
        while count['time'] <= t:
            one_step_evolution(count)
            if sum(count['pop_state'].values())==0:
                print('Population got extinct at time {}.'.format(count['time']))
                break
            elif 'stop' in count.keys():
                print('Evolution stopped at time {} with a total of {} individuals.'.format(count['time'],sum(count['pop_state'].values())))
                break
        if sum(count['pop_state'].values())==0:
            break
        if min(count['pop_state'].values()) < 0:
            print('One population state got negative')
            break
        #extract summerized statistics to safe in file
        stats_per_time.append(copy.deepcopy(extract_summerized_statistics(count)))
        print('Process: {}% at {}. Currently at time t={} there are {} individuals alive.'.format(t*100/model['t_max'],datetime.time(datetime.now()),t,sum(count['pop_state'].values())))
        
    t_end = time.clock()
    dt = t_end - t_start
      
    print('Evolution terminated successfully. At the final time t={} there were {} individuals alive. The simulation ran for {}.'.format(count['time'],sum(count['pop_state'].values()),display_time(dt)))
    print('The following Parameters have been used:')
    print(model)
    
    # Store list with locations after each round in a pickle file
    rel_path = "out_analysis/data/statistic_{}(K={},x_0=({}x{},{}),p^(-1)={},mating={}).pickle".format(
            model_name,
            model['K'],
            model['number_of_families'],
            int(model['K']*model["initial_rel_quantity_per_family"]),
            int(model['K']*model["initial_rel_quantity_of_disease"]),
            model['inv_mutation_rate'],
            model['mating_scheme']
            )
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "wb") as out_file:
        pickle.dump(stats_per_time, out_file)
