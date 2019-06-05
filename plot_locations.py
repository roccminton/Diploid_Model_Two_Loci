import json
import pickle
import matplotlib.pyplot as plt
import importlib
import os
import plotting_functions as pf
import seaborn as sns
 
def plot_something_over_time_one_axis(data, plotting_dict, t_0=0, t_max=None, step=1,safe=None):
    """Generic function to plot someting over time. The keys variable should be
    a list with first enty beeing the specification which plot one wants to realize
    and some additional arguments if necessary.
    
    Valid keys are *absolute*,*relative*,
    
    """
    
    fig, ax = plt.subplots()
    
    # Die obere und rechte Achse unsichtbar machen:
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    
    #get list indices handed in boundary times
    index_low, index_high = pf.set_time_index(t_0,t_max,stats_by_time)

    #Gives a sorted list of the timesteps
    time = pf.get_data(stats_by_time,'time',index_low,index_high,step)
    
    #set color map
    colors = sns.color_palette('hls', len(plotting_dict))
    ax.set_prop_cycle('color',colors)
        
    for specification, plotting_function in plotting_dict.items():
        plotting_function(data,specification,time,index_low,index_high,step,ax)
        

    # Die linke Diagrammachse auf den Bezugspunkt t_min legen:
    ax.spines['left'].set_position(('data',t_0))

    # Die untere Diagrammachse auf den Bezugspunkt '0' der y-Achse legen:
    ax.spines['bottom'].set_position(('data',0))
    
    #y-Achse immer bei 0 startend
    plt.ylim(bottom=0)
    
    plt.xscale('linear')
    
    #Labeling the Axes and adding titel
    plt.xlabel('Time')
    
    # Shrink current axis's height by factor of p on the top
    p = 0.025 * len(plotting_dict)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0,
                 box.width, box.height - box.height * p])

    #show legend
    ax.legend(loc='lower center', 
              bbox_to_anchor=(0.5,1),
              fontsize=10,
              edgecolor = 'None'
              )    
    
    if not safe==None:
        _safe_fig(safe,t_0,t_max)
    
    plt.show()

def plot_something_over_time_two_axis(data, plotting_dict1, plotting_dict2, t_0=0, t_max=None, step=1, safe=None, ylim=None):
    """Generic function to plot someting over time. The keys variable should be
    a list with first enty beeing the specification which plot one wants to realize
    and some additional arguments if necessary.
    
    Valid keys are *absolute*,*relative*,
    
    """
    fig, ax1 = plt.subplots()
        
    #get list indices handed in boundary times
    index_low, index_high = pf.set_time_index(t_0,t_max,stats_by_time)

    #Gives a sorted list of the timesteps
    time = pf.get_data(stats_by_time,'time',index_low,index_high,step)
    
    #set color map
    colors = sns.color_palette('BrBG', len(plotting_dict1))
    ax1.set_prop_cycle('color',colors)

        
    for specification, plotting_function in plotting_dict1.items():
            plotting_function(data,specification,time,index_low,index_high,step,ax1)
            
    #instantiate a second axes that shares the same x-axis        
    ax2 = ax1.twinx()
    
    #set color map
    colors = sns.color_palette('coolwarm', len(plotting_dict2))
    ax1.set_prop_cycle('color',colors)


    for specification, plotting_function in plotting_dict2.items():
            plotting_function(data,specification,time,index_low,index_high,step,ax2)
    
    plt.xscale('linear')
    
    #Die linke Diagrammachse auf den Bezugspunkt t_max legen:
    if t_max == None:
        t_max = model['t_max']

    for ax in ax1, ax2:
        #y-Achse immer bei 0 startend und nach oben begrenzt durch max Wert oder 1 wenn y_max kleiner ist
        ax.set_ylim(0)
        # Die obere und rechte Achse unsichtbar machen:
        ax.spines['top'].set_color('none')
        # Die linke Diagrammachse auf den Bezugspunkt t_min legen:
        ax.spines['left'].set_position(('data',t_0))
        # Die linke Diagrammachse auf den Bezugspunkt t_max legen:
        ax.spines['right'].set_position(('data',t_max))
        # Die untere Diagrammachse auf den Bezugspunkt '0' der y-Achse legen:
        ax.spines['bottom'].set_position(('data',0))
    
    #abs. Grenze für zweite Achse setzen
    if ylim != None:
        ax2.set_ylim(0,ylim)
        
    #show legend
    ax1.legend(
            loc='lower left', 
            bbox_to_anchor=(0,1),
            fontsize=9,
            edgecolor='None'
            )    
    
    #show legend
    ax2.legend(
            loc='lower right', 
            bbox_to_anchor=(1,1),
            fontsize=9,
            edgecolor='None'
            )    

    fig.tight_layout()
    
    #safe figure
    if not safe==None:
        _safe_fig(safe,t_0,t_max)

    plt.show()
   
def plot_something_over_time_three_axis(data, plotting_dict1, plotting_dict2, plotting_dict3, t_0=0, t_max=None, step=1, safe=None, ylim=None):
    """Generic function to plot someting over time. The keys variable should be
    a list with first enty beeing the specification which plot one wants to realize
    and some additional arguments if necessary. The first plotting_dict will be plotted
    on an invisible axis.
    

    """
    fig, ax1 = plt.subplots(figsize=(9.0,6.0))
        
    #get list indices handed in boundary times
    index_low, index_high = pf.set_time_index(t_0,t_max,stats_by_time)

    #Gives a sorted list of the timesteps
    time = pf.get_data(stats_by_time,'time',index_low,index_high,step)
    
    #set color map
    colors = sns.color_palette('BrBG', len(plotting_dict1))
    ax1.set_prop_cycle('color',colors)

        
    for specification, plotting_function in plotting_dict1.items():
            plotting_function(data,specification,time,index_low,index_high,step,ax1)
            
    #instantiate a second axes that shares the same x-axis        
    ax2 = ax1.twinx()
        
    #set color map
    colors = sns.color_palette('coolwarm', len(plotting_dict2))
    ax2.set_prop_cycle('color',colors)


    for specification, plotting_function in plotting_dict2.items():
            plotting_function(data,specification,time,index_low,index_high,step,ax2)
            
    #instantiate a third axes that shares the same x-axis        
    ax3 = ax2.twinx()
        
    #set color map
    colors = sns.color_palette('bright', len(plotting_dict3))
    ax3.set_prop_cycle('color',colors)


    for specification, plotting_function in plotting_dict3.items():
            plotting_function(data,specification,time,index_low,index_high,step,ax3)
    
    plt.xscale('linear')
    
    if t_max == None:
        t_max = model['t_max']

    
    for ax in ax1, ax2, ax3:
        #y-Achse immer bei 0 startend und nach oben begrenzt durch max Wert oder 1 wenn y_max kleiner ist
        ax.set_ylim(0)
        # Die obere und rechte Achse unsichtbar machen:
        ax.spines['top'].set_color('none')
        # Die linke Diagrammachse auf den Bezugspunkt t_min legen:
        ax.spines['left'].set_position(('data',t_0))
        # Die linke Diagrammachse auf den Bezugspunkt t_max legen:
        ax.spines['right'].set_position(('data',t_max))
        # Die untere Diagrammachse auf den Bezugspunkt '0' der y-Achse legen:
        ax.spines['bottom'].set_position(('data',0))
        
    #Achse 1 unsichtbar machen
    ax1.axes.get_yaxis().set_visible(False)
    
    #abs. Grenze für zweite Achse setzen
    if ylim != None:
        ax2.set_ylim(0,ylim)
        
    #show legend
    ax1.legend(
            loc='lower center', 
            bbox_to_anchor=(0.5,1),
            fontsize=9,
            edgecolor='None'
            )    
    
    #show legend
    ax2.legend(
            loc='lower left', 
            bbox_to_anchor=(0,1),
            fontsize=9,
            edgecolor='None'
            )    
    
    #show legend
    ax3.legend(
            loc='lower right', 
            bbox_to_anchor=(1,1),
            fontsize=9,
            edgecolor='None'
            )    

    fig.tight_layout()
    
    #safe figure
    if not safe==None:
        _safe_fig(safe,t_0,t_max)

    plt.show()

def _plot_rel_quantity(stats_by_time,specification,time,index_low,index_high,step,ax):
    """Plos the relative quantity of *healty*,homogeneous/heterogeneous individuals in all homogeneous/heterogeneous
    individuals and the relative quantity of individuals carrying the illnes.
    
    """    
    #safe specification tuple in single variables
    h_nom, m_nom, h_denom, m_denom = specification
        
    #setup label_nominator and label_denominator list for the plotting function
    label_nominator = _setup_label_list(h_nom,m_nom,stats_by_time,index_high)
    label_denominator = _setup_label_list(h_denom,m_denom,stats_by_time,index_high)
        
    return ax.plot(
            time,
            pf.get_rel_quantities(stats_by_time,label_nominator,label_denominator,index_low, index_high, step),
            label = _get_label_from_specification(specification)
            )
        
def _plot_abs_quantity(stats_by_time,specification,time,index_low,index_high,step,ax):
    """Plos the relative quantity of *healty*,homogeneous/heterogeneous individuals in all homogeneous/heterogeneous
    individuals and the relative quantity of individuals carrying the illnes.
    
    """
    #safe specification tuple in single variables
    h,m = specification
            
    #setup label list
    label_list = _setup_label_list(h,m,stats_by_time,index_high)
    
    if len(label_list) == 0:
        raise ValueError('The lable {} is not represented in the current population'.format(_get_label_from_specification(specification)))
    else:
        return ax.plot(
                time,
                pf.get_abs_quantities(stats_by_time,label_list,index_low, index_high, step),
                label = _get_label_from_specification(specification)
                )

def _plot_data(stats_by_time,specification,time,index_low,index_high,step,ax):
    """Plots the data over the time interval [index_low,index_high] with given stepsize.
    Valid data types are ['time','trait','population_size','mutation_counter','mutation_load'].
    The specification must consist of a tuple of the following form
        (data_type, label, color)
    
    """
    data_type, label, color, marker = specification
    
    #Plot data line over time in given color
    return ax.plot(
            time,
            pf.get_data(stats_by_time,data_type,index_low,index_high,step),
            label=label,
            color=color,
            marker=marker,
            markevery=int(model['t_max']*0.1)
            )
    
def _plot_data_over_pop_size(stats_by_time,specification,time,index_low,index_high,step,ax):
    """Plots the data over the time interval [index_low,index_high] with given stepsize.
    Valid data types are ['time','trait','population_size','mutation_counter','mutation_load'].
    The specification must consist of a tuple of the following form
        (data_type, label, color)
    
    """
    data_type, label, color, marker = specification
    
    #Plot data line over time in given color
    return ax.plot(
            time,
            pf.get_data_over_pop_size(stats_by_time,data_type,index_low,index_high,step),
            label=label,
            color=color,
            marker=marker,
            markevery=int(model['t_max']*0.1)
            )

def _plot_rel_entropy(stats_by_time,specification,time,index_low,index_high,step,ax):
    """Plots the relative entropy over the time interval [index_low,index_high] with given stepsize.
    """
    data_type, label, color, marker = specification
    
    #Plot data line over time in given color
    return ax.plot(
            time,
            pf.get_rel_entropy_of_family_sizes(stats_by_time,index_low,index_high,step),
            label=label,
            color=color,
            marker=marker,
            markevery=int(model['t_max']*0.1)
            )


def _setup_label_list(h,m,stats_by_time,index_high):
    """Input must be of the form (h,m), where h in {'homo','hetero','total'} and
    m in {0,1,'total'}
    
    """
    #save pop_state keys in list
    pop_state_keys = list(stats_by_time[index_high-1]['pop_state'].keys())
        
    #if the homo/hetero specification is 'total'
    if h == 'total':
        #if the mutation/no mutation specification is 'total'
        if m == 'total':
            return pop_state_keys
        elif m == '>0':
            return list(filter(lambda key: key[1]>0, pop_state_keys))
        elif type(m) == int and m >= 0:
            return list(filter(lambda key: key[1]==m, pop_state_keys))
        else:
            raise KeyError('The mutation specification {} is invalide.'.format(m))
    elif h in {'homo','hetero'}:
        #if the mutation/no mutation specification is 'total'
        if m == 'total':
            return list(filter(lambda key: key[0]==h, pop_state_keys))
        elif m == '>0':
            return list(filter(lambda key: key[0]==h and key[1]>0, pop_state_keys))
        elif type(m) == int and m >= 0:
            return list(filter(lambda key: key[0]==h and key[1]==m, pop_state_keys))
        else:
            raise KeyError('The mutation specification {} is invalide.'.format(m))
    else:
        raise KeyError('The homogeneous/heterogeneous specification {} is invalide.'.format(h))
        
def _get_label_from_specification(specification):
    """Returns the lable for the plot based on the specification variable.
    
    """
    if type(specification) == tuple:
        #the case of absolute plots
        if len(specification) == 2:
            if specification[0] == 'total':
                if specification[1] == 'total':
                    return '# total population'
                elif specification[1] == '>0':
                    return '# of carrier of the disease'
                elif type(specification[1]) == int: 
                    if specification[1] > 0:
                        return '# carrying the mutation {}'.format(_get_word_from_integer(specification[1]))
                    elif specification[1] == 0:
                        return '# not carrying the mutation'
                else:
                    raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
            elif specification[0] in ['homo','hetero']:
                if specification[1] == 'total':
                    return '# {}geneous'.format(specification[0])
                elif specification[1] == '>0':
                    return '# {}geneous carrier of the disease'.format(specification[0])
                elif type(specification[1]) == int: 
                    if specification[1] > 0:
                        return '# of {}geneous carrying the mutation {}'.format(specification[0],_get_word_from_integer(specification[1]))
                    elif specification[1] == 0:
                        return '# of {}geneous not carrying the mutation'.format(specification[0])
                    else:
                        raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                else:
                    raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
        #the case of relative plots
        elif len(specification) == 4:
            if specification[0] == 'total':
                if specification[1] == '>0':
                    if specification[2] == specification[3] == 'total':
                        return '% of carrier of disease'
                    else:
                        raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                elif type(specification[1]) == int:
                    if specification[1] > 0:
                        if specification[3]  == '>0':
                            return '% of carrier of disease carrying the mutation {}'.format(_get_word_from_integer(specification[1]))
                        elif specification[3] == 'total':
                            return '% carrying the mutation {}'.format(_get_word_from_integer(specification[1]))
                        else:
                            raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                    elif specification[1] == 0:
                        return '% not carrying the mutation'
                    else:
                        raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                else:
                    raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
            elif specification[0] in ['homo','hetero']:
                if specification[1] == 'total':
                    if specification[2] == specification[3] == 'total':
                        return '% of {}geneous'.format(specification[0])
                    else:
                        raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                elif specification[1] == '>0':
                    if specification[2] == 'total':
                        if specification[3] == 'total':
                            return '% of {}geneous carrier of disease'.format(specification[0])
                        elif specification[3] == specification[1]:
                            return '% of {}geneous carrier of disease in all carrier'.format(specification[0])
                        else:
                            raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                    elif specification[2] == specification[0] and specification[3] == 'total':
                        return '% of {}geneous carrier of disease in all {}geneous.'.format(specification[0],specification[0])
                    else:
                        raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))                        
                elif type(specification[1]) == int:
                    if specification[1] > 0:
                        if specification[2] == 'total':
                            if specification[3] == 'total':
                                return '% of {}geneous carrying the mutation {}'.format(specification[0],_get_word_from_integer(specification[1]))
                            elif specification[3] == '>0':
                                return '% of {}geneous carrying the mutation {} in all carrier of disease'.format(specification[0],_get_word_from_integer(specification[1]))
                            elif specification[3] == specification[1]:
                                return '% of {}genous in all carrying the mutation {}'.format(specification[0],_get_word_from_integer(specification[1]))
                            else:
                                raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                        elif specification[2] == specification[0]:
                            if specification[3] == 'total':
                                return '% of carrying the mutation {} in all {}geneous'.format(specification[0])
                            elif specification[3] == '>0':
                                return '% of {}geneous carrying the mutation {} in all carrier of disease'.format(specification[0],_get_word_from_integer(specification[1]))
                            else:
                                raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                        else:
                            raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                    elif specification[1] == 0:
                        if specification[2] == 'total':
                            if specification[3] == 'total':
                                return '% of {}geneous not carrying the mutation'.format(specification[0])
                            elif specification[3] == specification[1]:
                                return '% of not carrying the mutation in all {}geneous'.format(specification[0])
                            else:
                                raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                        elif specification[2] == specification[0]:
                            if specification [3] == 'total':
                                return '% of not carrying the mutation in all {}geneous'.format(specification[0])
                            else:
                                raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                        else:
                            raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                    else:
                        raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
                else:
                    raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
            else:
                raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
        else:
            raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))
    else:
        raise ValueError('Cannot convert the specification {} to a lable.'.format(specification))

def _get_word_from_integer(n):
    """Returns a string with the word for the frequency of the integer.
    
    """
    #convert to integer
    try:
        n = int(n)
    except TypeError or ValueError:
        raise ValueError('Cannot convert {} into word'.format(n))
    #get word
    if n not in [1,2]:
        return '{}-times'.format(n)
    elif n == 1:
        return 'once'
    elif n == 2:
        return 'twice'
    else:
        raise ValueError('Cannot convert {} into word'.format(n))
        
def _safe_fig(title,t_0,t_max):
    """Safes the figure with the respective title.
    
    """
    if model['mating_scheme'][0] == 'dynasty':
        mating_scheme = 'dynasty(nmr={})'.format(model['mating_scheme'][1])
    elif model['mating_scheme'][0] == 'random':
        mating_scheme = 'random'
    elif model['mating_scheme'][0] == 'western':
        mating_scheme = 'western(mr={})'.format(model['mating_scheme'][1])
    else:
        raise ValueError('Mating scheme *{}* in _safe_fig unknown.'.format(model['mating_scheme'][0]))
            
    initial_state = '({}x{},{})'.format(model['number_of_families'],int(model['K']*model['initial_rel_quantity_per_family']),int(model['K']*model['initial_rel_mutation_load']))
    path = "{}_{}(K={},x_0={},p^(-1)={}".format(title,mating_scheme,model['K'],initial_state,model['inv_mutation_rate'])
    
    if t_0>0:
        path += ',t_0={}'.format(t_0)
    if t_max != None:
        path += ',t_max={}'.format(t_max)
        
    path += ').pdf'
    #safe image 
    plt.savefig(os.path.join(r'/home/larocca/Dokumente/Simulationen/Diploid_Model_Two_Loci/out_analysis/figures',
                                 path))
        
if __name__ == "__main__":
    #Remember to manually set the model name also in plotting_functions.py
    model_name = 'Families'
    
    #Absolute dir the script is in. Needs to be changed manually if the file is moveed
    script_dir = r'/home/larocca/Dokumente/Simulationen/Diploid_Model_Two_Loci'
    
    #load model specific variables such as carrying capacity, initial trait and marker, etc.
    rel_path = "model_specifications/{}.json".format(model_name)
    abs_file_path = os.path.join(script_dir, rel_path)
    model = json.load(open(abs_file_path), encoding="utf-8")
    
#    #set model parameters manually
#    model['K'] = 100000
#    model['initial_rel_mutation_load'] = 0.0517
#    model['initial_rel_quantity_per_family'] = 0.0002
#    model['inv_mutation_rate'] = 1200
    model['mating_scheme'] = ['random']
    model['t_max'] = 10000
    model['number_of_loci'] = 1000
    
    # Load locations from pickle file
    rel_path = "out_analysis/data/statistic_{}(K={},nloci={},x_0=({}x{},{}),p^(-1)={},mating={}).pickle".format(
            model_name,
            model['K'],
            model['number_of_loci'],
            model['number_of_families'],
            int(model['K']*model["initial_rel_quantity_per_family"]),
            int(model['K']*model["initial_rel_mutation_load"]),
            model['inv_mutation_rate'],
            model['mating_scheme']
            )
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "rb") as in_file:
        stats_by_time = pickle.load(in_file)

    #import model specific functions such as birthrate, mutation operator and probabilities, etc.
    spc = importlib.import_module('model_specifications.{}'.format(model_name))
    
    plot_dict_abs = {
            ('total','total'):_plot_abs_quantity,
#            ('total','>0'):_plot_abs_quantity,
#            ('total',3):_plot_abs_quantity,
#            ('total',2):_plot_abs_quantity,
#            ('total',1):_plot_abs_quantity,
#            ('total',0):_plot_abs_quantity,
#            ('homo', 'total'):_plot_abs_quantity,
#            ('homo', 1):_plot_abs_quantity,
#            ('homo', 0):_plot_abs_quantity,
#            ('hetero','total'):_plot_abs_quantity,
#            ('hetero',1):_plot_abs_quantity,
#            ('hetero',0):_plot_abs_quantity,
            }
    
    plot_dict_rel = {
#            ('total',4,'total','total'):_plot_rel_quantity,            
#            ('total',3,'total','total'):_plot_rel_quantity,
#            ('total',2,'total','total'):_plot_rel_quantity,
#            ('total',1,'total','total'):_plot_rel_quantity,
            ('total','>0','total','total'):_plot_rel_quantity,
            ('homo','total','total','total'):_plot_rel_quantity,
#            ('homo',1,'total','total'):_plot_rel_quantity,
#            ('homo',1,'homo','total'):_plot_rel_quantity,
#            ('homo',1,'total',1):_plot_rel_quantity,
#            ('homo',0,'total','total'):_plot_rel_quantity,
#            ('homo',0,'homo','total'):_plot_rel_quantity,
#            ('homo',0,'total',0):_plot_rel_quantity,
#            ('hetero','total','total','total'):_plot_rel_quantity,
#            ('hetero',1,'total','total'):_plot_rel_quantity,
#            ('hetero',1,'hetero','total'):_plot_rel_quantity,
#            ('hetero',1,'total',1):_plot_rel_quantity,
#            ('hetero',0,'total','total'):_plot_rel_quantity,
#            ('hetero',0,'hetero','total'):_plot_rel_quantity,
#            ('hetero',0,'total',0):_plot_rel_quantity            
        }
    
    plot_dict_data = {
#            ('mutation_load','Mutation load','r'):_plot_data,
            ('number_false', 'Number of ill individual','peru',None):_plot_data,
#            ('number_of_types','Number of different family-types','g',None): _plot_data,
#            ('mutation_counter', 'Total umber of mutation','g'):_plot_data,
#             ('active_families','# families alive','lightgrey'):_plot_data,
            }
    
    plot_dict_data_rel = {
            ('mutation_load','Relative mutation load','r',None):_plot_data_over_pop_size,
            ('bla','Relative Entropy of family sizes','b',None):_plot_rel_entropy,
            }

#    title = 'rel_entropy'
#    title='%_indv_mut_load'
#    title='mutation_load'
    title = 'ill_individal'
    
    plot_something_over_time_three_axis(stats_by_time,plot_dict_abs,plot_dict_data,plot_dict_data_rel,safe=title,t_max=5000)
#    plot_something_over_time_two_axis(stats_by_time,plot_dict_data,plot_dict_data_rel)
#    plot_something_over_time_one_axis(stats_by_time,plot_dict_data_rel)
#    plot_something_over_time_one_axis(stats_by_time,plot_dict_abs,t_max=300)    