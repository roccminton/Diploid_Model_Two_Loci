import json
import pickle
import matplotlib.pyplot as plt
import importlib
import os
import plotting_functions as pf
 
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

def plot_something_over_time_two_axis(data, plotting_dict1, plotting_dict2, t_0=0, t_max=None, step=1, safe=None):
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
        
    for specification, plotting_function in plotting_dict1.items():
            plotting_function(data,specification,time,index_low,index_high,step,ax1)
            
    #instantiate a second axes that shares the same x-axis        
    ax2 = ax1.twinx()
    
    # Die obere und rechte Achse unsichtbar machen:
    ax2.spines['top'].set_color('none')

    for specification, plotting_function in plotting_dict2.items():
            plotting_function(data,specification,time,index_low,index_high,step,ax2)
    
    plt.xscale('linear')
    
    #y-Achse immer bei 0 startend und nach oben begrenzt durch max Wert oder 1 wenn y_max kleiner ist
    ax1.set_ylim(0,max(ax1.get_ylim()[1],1))
    ax2.set_ylim(0,max(ax2.get_ylim()[1],1))
    
    # Die obere und rechte Achse unsichtbar machen:
    ax1.spines['top'].set_color('none')
    ax2.spines['top'].set_color('none')
    # Die linke Diagrammachse auf den Bezugspunkt t_min legen:
    ax1.spines['left'].set_position(('data',t_0))
    ax2.spines['left'].set_position(('data',t_0))
    # Die linke Diagrammachse auf den Bezugspunkt t_max legen:
    if t_max == None:
        t_max = model['t_max']
    ax1.spines['right'].set_position(('data',t_max))
    ax2.spines['right'].set_position(('data',t_max))
    # Die untere Diagrammachse auf den Bezugspunkt '0' der y-Achse legen:
    ax1.spines['bottom'].set_position(('data',0))
    ax2.spines['bottom'].set_position(('data',0))
        
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
   

def _plot_rel_quantity(stats_by_time,specification,time,index_low,index_high,step,ax):
    """Plos the relative quantity of *healty*,homogeneous/heterogeneous individuals in all homogeneous/heterogeneous
    individuals and the relative quantity of individuals carrying the illnes.
    
    """    
    #safe specification tuple in single variables
    h_nom, m_nom, h_denom, m_denom = specification

    #set color and label dictionary
    #set offset, onoff sequence for dotted lines
    dotted = [1,1]      
    dashed = [5,5]
    dashdot = [1,5,3,5]    #i.e. 1pt line, 5pt break, 3pt line, 5pt break
    
    color_and_label = {
            ('total',1,'total','total'):['blue','% of (0,1)',dotted],
            ('total',0,'total','total'):['orange','% of (0,0)',dotted],
            ('homo','total','total','total'):['red','% of homogeneous',dotted],
            ('homo',1,'total','total'):['purple','total % of [homogeneous, (0,1)]',dotted],
            ('homo',1,'homo','total'):['purple','% of [homogenous, (0,1)] in all homogeneous',dotted],
            ('homo',1,'total',1):['purple','% of [homogenous, (0,1)] in all (0,1)s',dotted],
            ('homo',0,'total','total'):['darkorange','total % of [homogeneous, (0,0)]',dotted],
            ('homo',0,'homo','total'):['darkorange','% of [homogenous, (0,0)] in all homogeneous',dotted],
            ('homo',0,'total',0):['darkorange','% of [homogenous, (0,0)] in all (0,0)',dotted],
            ('hetero','total','total','total'):['green','% of heterogeneous',dotted],
            ('hetero',1,'total','total'):['c','total % of [heterogeneous, (0,1)]',dotted],
            ('hetero',1,'hetero','total'):['c','% of [heterogenous, (0,1)] in all heterogeneous',dotted],
            ('hetero',1,'total',1):['c','% of [heterogenous, (0,1)] in all (0,1)s',dotted],
            ('hetero',0,'total','total'):['gold','total % of [heterogeneous, (0,0)]',dotted],
            ('hetero',0,'hetero','total'):['gold','% of [heterogenous, (0,0)] in all heterogeneous',dotted],
            ('hetero',0,'total',0):['gold','% of [heterogenous, (0,0)] in all (0,0)',dotted]            
        }
        
    #setup label_nominator and label_denominator list for the plotting function
    label_nominator = _setup_label_list(h_nom,m_nom)
    label_denominator = _setup_label_list(h_denom,m_denom)
        
    return ax.plot(
            time,
            pf.get_rel_quantities(stats_by_time,label_nominator,label_denominator,index_low, index_high, step),
            color = color_and_label[specification][0],
            label = color_and_label[specification][1],
            dashes = color_and_label[specification][2],
            )
        
def _plot_abs_quantity(stats_by_time,specification,time,index_low,index_high,step,ax):
    """Plos the relative quantity of *healty*,homogeneous/heterogeneous individuals in all homogeneous/heterogeneous
    individuals and the relative quantity of individuals carrying the illnes.
    
    """
    #safe specification tuple in single variables
    h,m = specification
    
    #set color and label dictionary
    color_and_label = {
            'total':{
                    'total':['grey','# Total population'],
                    1:['blue','# (0,1)'],
                    0:['orange','# (0,0)']
                    },
            'homo': {
                    'total':['red','# homogeneous'],
                    1:['purple','# [homogeneous, (0,1)]'],
                    0:['darkorange','# [homogeneous, (0,0)]']
                    },
            'hetero': {
                    'total':['green','# heterogeneous'],
                    1:['c','# [heterogeneous, (0,1)]'],
                    0:['gold','# [heterogeneous, (0,0)]']
                    }
        }
        
    #setup label list
    label_list = _setup_label_list(h,m)
    
    return ax.plot(
            time,
            pf.get_abs_quantities(stats_by_time,label_list,index_low, index_high, step),
            color = color_and_label[h][m][0],
            label = color_and_label[h][m][1]
            )

def _plot_number_of_active_families(stats_by_time,specification,time,index_low,index_high,step,ax):
    """Plots the number of active families over the time interval
    [index_low,index_high] with given stepsize.
    
    """
    #Plot data line over time in given color
    return ax.plot(
            time,
            pf.get_number_of_active_families(stats_by_time,index_low,index_high,step),
            label='# families alive',
            color='lightgrey',
            marker='.',
            markevery=int(model['t_max']*0.1)
            )

def _setup_label_list(h,m):
    """Input must be of the form (h,m), where h in {'homo','hetero','total'} and
    m in {0,1,'total'}
    
    """    
    #if the homo/hetero specification is 'total'
    if h == 'total':
        #if the mutation/no mutation specification is 'total'
        if m == 'total':
            return [('homo',0),('homo',1),('hetero',0),('hetero',1)]
        elif m in {0,1}:
            return [('homo',m),('hetero',m)]
        else:
            raise KeyError('The mutation specification {} is invalide.'.format(m))
    elif h in {'homo','hetero'}:
        #if the mutation/no mutation specification is 'total'
        if m == 'total':
            return [(h,0),(h,1)]
        elif m in {0,1}:
            return [(h,m)]
        else:
            raise KeyError('The mutation specification {} is invalide.'.format(m))
    else:
        raise KeyError('The homogeneous/heterogeneous specification {} is invalide.'.format(h))
        
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
            
    initial_state = '({}x{},{})'.format(model['number_of_families'],int(model['K']*model['initial_rel_quantity_per_family']),int(model['K']*model['initial_rel_quantity_of_disease']))
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
    model['mating_scheme'] = ['dynasty', 0.1]
    
    # Load locations from pickle file
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
    with open(abs_file_path, "rb") as in_file:
        stats_by_time = pickle.load(in_file)

    #import model specific functions such as birthrate, mutation operator and probabilities, etc.
    spc = importlib.import_module('model_specifications.{}'.format(model_name))
    
    plot_dict_abs = {
#            ('total','total'):_plot_abs_quantity,
#            ('total',1):_plot_abs_quantity,
#            ('total',0):_plot_abs_quantity,
            ('homo', 'total'):_plot_abs_quantity,
#            ('homo', 1):_plot_abs_quantity,
#            ('homo', 0):_plot_abs_quantity,
            ('hetero','total'):_plot_abs_quantity,
#            ('hetero',1):_plot_abs_quantity,
#            ('hetero',0):_plot_abs_quantity,
            }
    
    plot_dict_rel = {
            ('total',1,'total','total'):_plot_rel_quantity,
#            ('total',0,'total','total'):_plot_rel_quantity,
#            ('homo','total','total','total'):_plot_rel_quantity,
#            ('homo',1,'total','total'):_plot_rel_quantity,
            ('homo',1,'homo','total'):_plot_rel_quantity,
#            ('homo',1,'total',1):_plot_rel_quantity,
#            ('homo',0,'total','total'):_plot_rel_quantity,
#            ('homo',0,'homo','total'):_plot_rel_quantity,
#            ('homo',0,'total',0):_plot_rel_quantity,
#            ('hetero','total','total','total'):_plot_rel_quantity,
#            ('hetero',1,'total','total'):_plot_rel_quantity,
            ('hetero',1,'hetero','total'):_plot_rel_quantity,
#            ('hetero',1,'total',1):_plot_rel_quantity,
#            ('hetero',0,'total','total'):_plot_rel_quantity,
#            ('hetero',0,'hetero','total'):_plot_rel_quantity,
#            ('hetero',0,'total',0):_plot_rel_quantity            
        }

    family_dic = {'something':_plot_number_of_active_families}
    
    title='%of(0,1)'
    
    plot_something_over_time_two_axis(stats_by_time,plot_dict_abs,family_dic,safe=None)
#    plot_something_over_time_one_axis(stats_by_time,plot_dict_rel)
#    plot_something_over_time_one_axis(stats_by_time,plot_dict_abs,t_max=300)    