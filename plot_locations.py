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
        _safe_fig(safe,t_0,t_max,step)

    plt.show()
   
def plot_something_over_time_three_axis(data, plotting_dict1, 
                                        plotting_dict2, 
                                        plotting_dict3, 
                                        t_0=0, 
                                        t_max=None, 
                                        step=1, 
                                        safe=None, 
                                        ylim=None,
                                        family_sizes = None):
    """Generic function to plot someting over time. The keys variable should be
    a list with first enty beeing the specification which plot one wants to realize
    and some additional arguments if necessary. The first plotting_dict will be plotted
    on an invisible axis.
    

    """
    
    a_mut = 1
    a_ill = 1
    
    #size_factor
    s = 1
    
    fig, ax1 = plt.subplots(figsize=(s*9.0,s*6.0))
        
    #get list indices handed in boundary times
    index_low, index_high = pf.set_time_index(t_0,t_max,stats_by_time)

    #Gives a sorted list of the timesteps
#    time = pf.get_data(stats_by_time,'time',index_low,index_high,step)
    time = list(range(index_low, index_high, step))
    
    #set x-lable to mating scheme
#    x_label = 'Mating scheme: {}'.format(_get_name_from_mating_scheme())
    if family_sizes != None:
        x_label = '\n average family size: {}'.format(family_sizes[1])
    ax1.set_xlabel(x_label)
    
    #set color map
    colors = sns.color_palette('BrBG', len(plotting_dict1))
    ax1.set_prop_cycle('color',colors)
#    ax1.set_prop_cycle(color=['k'],alpha=[.5],linewidth=[s])
        
    for specification, plotting_function in plotting_dict1.items():
            plotting_function(data,specification,time,index_low,index_high,step,ax1)
            
    #instantiate a second axes that shares the same x-axis        
    ax2 = ax1.twinx()
        
    #set color map
    colors = sns.color_palette('coolwarm', len(plotting_dict2))
    ax2.set_prop_cycle('color',colors)
#    ax2.set_prop_cycle(alpha=[a_ill],linewidth=[s])


    for specification, plotting_function in plotting_dict2.items():
            plotting_function(data,specification,time,index_low,index_high,step,ax2)
            
    #instantiate a third axes that shares the same x-axis        
    ax3 = ax2.twinx()
        
    #set color map
    colors = sns.color_palette('bright', len(plotting_dict3))
    ax3.set_prop_cycle('color',colors)
#    ax3.set_prop_cycle(alpha=[a_mut],linewidth=[s])

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
        # untere Achse keine ticks
#        ax.set_xticks([])
    
#        
#    ax1.set_xlabel('time',fontsize=11*s)
#    ax1.xaxis.set_label_coords(0.95,-0.015)
        
#    Achse 1 unsichtbar machen
    ax1.axes.get_yaxis().set_visible(False)
    
    #Total Population Achsenbeschriftung
#    a_total = 0
#    ax1.set_yticks([500,10000])
#    ax1.set_yticklabels(['500','10.000'], alpha=a_total, fontsize=11*s)
#    ax1.tick_params(axis='y',length=0, pad=10)

    #add vertical line for 'now'
#    a_now = 0 #alpha wert für an/aus
#    t = 15500
#    ax1.axvline(x=t, color='r',alpha=a_now, linewidth = 2*s)
#    ax1.set_xticks([t-(10000+2020),t-2020,t+(10000-2020),t+(20000-2020),t+(30000-2020)])
#    ax1.set_xticklabels(['-10.000','0','10.000','20.000','30.000'], color='k',alpha=1, fontsize=11*s)
#    ax1.tick_params(axis='x',length=9, pad=10)
#    ax1.text(t-850,-300,"now",fontsize=11*s, color='r',alpha=a_now)
    
    #rechte Achse
#    color_ill = list(plotting_dict2.keys())[0][2]
    ax2.set_ylim(0,0.053)
#    ax2.set_yticks([0.01,0.02,0.03,0.04,0.05])
#    ax2.set_yticklabels(['1%','2%','3%','4%','5%'], alpha=a_ill, fontsize=13*s, color=color_ill)
#    ax2.tick_params(axis='y',length=0, pad=10)
#    ax1.spines['left'].set_color('none')
#    ax3.spines['left'].set_color('none')
#    ax2.spines['left'].set_color(color_ill)
#    ax2.spines['left'].set_linewidth(3)
    
    #rechte Achse (für Differenzenplot)
#    ax2.set_ylim(-0.0325,0.0325)
#    ax2.set_yticks([-0.03,-0.02,-0.01,0.01,0.02,0.03])
#    ax2.set_yticklabels(['-3%','-2%','-1%','1%','2%','3%'], alpha=a_ill, fontsize=11*s, color=color_ill)
#    ax2.tick_params(axis='y',length=0, pad=10)

    #rechte Achse (für Quotient)
#    ax2.set_ylim(0,3.1)
#    ax2.set_yticks([1,2,3])
#    ax2.set_yticklabels(['1x','2x','3x'], alpha=a_ill, fontsize=13*s, color=color_ill)
#    ax2.tick_params(axis='y',length=0, pad=10)
    
    #linke Achse
#    color_mut = list(plotting_dict3.keys())[0][2]
    ax3.set_ylim(0,8.53)
#    ax3.set_yticks([1,2,3,4,5,6,7,8])
#    ax3.set_yticklabels(['1','2','3','4','5','6','7','8'], alpha=a_mut, fontsize=13*s, color=color_mut)
#    ax3.tick_params(axis='y',length=0, pad=10)
#    if a_mut != 0:
#        ax3.spines['right'].set_color(color_mut)
#        ax3.spines['right'].set_linewidth(3)
    
    #abs. Grenze für zweite Achse setzen
    if ylim != None:
        ax2.set_ylim(0,ylim)
        
    #show legend
#    ax1.legend(
#            loc='lower center', 
#            bbox_to_anchor=(0.5,1),
#            fontsize=s*9,
#            edgecolor='None'
#            )    
#    
#    if a_ill != 0:
#        #show legend
#        ax2.legend(
#                loc='lower left', 
#                bbox_to_anchor=(0,1),
#                fontsize=s*9,
#                edgecolor='None'
#                )    
#    
#    if a_mut != 0:
#        #show legend
#        ax3.legend(
#                loc='lower right', 
#                bbox_to_anchor=(1,1),
#                fontsize=s*9,
#                edgecolor='None'
#                )    

    fig.tight_layout()
    
    #safe figure
#    if not safe==None:
#        _safe_fig(safe,t_0,t_max,step,family_sizes)
    
    path = "Consang_family_size{}(same_ax).pdf".format(family_sizes)
    
    plt.savefig(os.path.join(r'/home/larocca/Dokumente/Simulationen/Diploid_Model_Two_Loci/out_analysis/Paper_Krawitz',
                                 path), transparent=True)

    plt.show()
    
def _add_average_line(ax,data,data_min,data_max,xmin,xmax,color='k'):
    """Adds the average of the data list from data_min to data_max as a horizontal line
    from xmin to xmax in the handed in color to the ax plot.
    
    """
    average = _get_average(data,data_min,data_max)
    ax.text(1/20,9/10,'{}'.format(round(average,2)),transform=ax.transAxes)
    ax.hlines(average,xmin,xmax,color)

def _get_average(data,data_min,data_max):
    """Calculates the average of a list over an index range.
    
    """
    
    return sum(data[data_min:data_max])/len(data[data_min:data_max])

def _get_averages_different_parameters(parameters_K_tmax,parameters_n_loci,data_type):
    """Retruns a dictionary with (x,y):z entries where z is the averaged parameter.
    
    """
    total = len(parameters_K_tmax)*len(parameters_n_loci)
    run = 0
    
    averages = {}
    for K, t_max in parameters_K_tmax:
        for n_loci in parameters_n_loci:
            #    #set model parameters manually
            model['K_list'] = [K]
            model['t_max'] = t_max
            model['number_of_loci'] = n_loci

            #set the relativ path
            rel_path = "out_analysis/data/statistic_{}(K={},nloci={},mating={},n_runs={}".format(
                    model_name,
                    model['K_list'],
                    model['number_of_loci'],
                    model['mating_scheme'],
                    model['num_runs']
                    )
            if model['mating_scheme'][0] == 'consang':
                rel_path += ",family_sizes=[100]).pickle"
            else:
                rel_path += ").pickle"
            # Load locations from pickle file
            abs_file_path = os.path.join(script_dir, rel_path)
            try:
                with open(abs_file_path, "rb") as in_file:
                    stats_by_time = pickle.load(in_file)
                    
                averages[(K,n_loci)] = _get_average(pf.get_data_over_pop_size(stats_by_time,data_type),int(t_max/2),int(t_max))
            except FileNotFoundError:
                averages[(K,n_loci)] = None
            
            #Print status to console
            run += 1
            print('Loading data: {}% ...please wait.'.format(int(100*(run/total))))
            
    print('Loading complete.')
            
    return averages

def plot_averages(averages_dict,x_axis,data_type,safe=False):
    """Plots the generated averages on the y axis and the different equ. pop_sizes on the 
    x axis. Different n_loci are displayed as different colors if the x_axis lable is set to pop_size.
    Instead if x_axis is n_loci it is the other way around.
    
    """
    if x_axis == 'pop_size':
        index = 0
        color_palette = 'GnBu_d'
        axis_title = "Equilibrium Population Size"
        legend_title = 'Number of Genes'
        tick_rotation = 30
    elif x_axis == 'n_loci':
        index = 1
        color_palette = 'OrRd_d'
        axis_title = "Number of Genes"
        legend_title = 'Equilibrium Population Size'
        tick_rotation = 30
    else:
        raise ValueError('Handed in x Lable is unknown: {}'.format(x_axis))
    
    #setup the figure
    fig, ax = plt.subplots(figsize=(9.0,6.0))
    # Die obere und rechte Achse unsichtbar machen:
    ax.spines['top'].set_color('none')
    # Die untere Diagrammachse auf den Bezugspunkt '0' der y-Achse legen:
    ax.spines['bottom'].set_position(('data',0))
    #Titel setzten
    ax.set_xlabel(axis_title)
    if data_type == 'number_false':
        ax.set_ylabel('Equilibrium % of ill individual')
    elif data_type == 'mutation_load':
        ax.set_ylabel('Equilibrium Relative Mutation Load')
    else:
        raise ValueError('Add y-Axis lable for the new data_type *{}* before plotting.'.format(data_type))
    
    #get the set of all parameters
    parameter_set = set()
    for K_and_n_loci in averages_dict.keys():
        parameter_set.add(K_and_n_loci[1-index])
        
    #set color map
    colors = sns.color_palette(color_palette, len(parameter_set))
    ax.set_prop_cycle('color',colors)
    
    for parameter in sorted(list(parameter_set)):
        
        #build up the lists for the plot
        xdata = []
        ydata = []
        
        for K_and_n_loci, average in averages_dict.items():
            if K_and_n_loci[1-index] == parameter:
                xdata.append(K_and_n_loci[index])
                ydata.append(average)
                    
        #add the plot to the figure
        ax.plot(
                xdata,
                ydata,
                label=parameter
                )
        
    #set the x-axis ticks
    #don't show the tick at pop_size 1000 because it is too close betwee 500-2000
    #and would be covered by the sorrounding ticks
    if x_axis == 'pop_size':
        xdata.remove(1000)
    plt.xticks(ticks=xdata,rotation=tick_rotation)

    #show legend
    ax.legend(
            loc='upper left', 
            bbox_to_anchor=(0,1),
            fontsize=9,
            edgecolor='k',
            title = legend_title
            )        
    
    if safe:
        path = "Eq_{}_vs_{}_{}.pdf".format(data_type,x_axis,_get_name_from_mating_scheme())
        #safe image 
        plt.savefig(os.path.join(r'/home/larocca/Dokumente/Simulationen/Diploid_Model_Two_Loci/out_analysis/figures',
                                     path))

        
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
    
    data = pf.get_rel_quantities(stats_by_time,label_nominator,label_denominator,index_low, index_high, step)
            
    return ax.plot(
            time,
            data,
            label = _get_label_from_specification(specification)
            )
        
def _plot_abs_quantity(stats_by_time,specification,time,index_low,index_high,step,ax,average_line=False):
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
    data_type, label, color, marker, average_line = specification
    
    data = pf.get_data(stats_by_time,data_type,index_low,index_high,step)
    
    if average_line:
        _add_average_line(
                ax,
                data,
                int(len(data)/2),
                int(len(data)),
                time[0],
                time[-1])               

    
    #Plot data line over time in given color
    return ax.plot(
            time,
            data,
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
    data_type, label, color, marker, average_line = specification
    
    data = pf.get_data_over_pop_size(stats_by_time,data_type,index_low,index_high,step)
    
    if average_line:
        _add_average_line(
                ax,
                data,
                int(len(data)/2),
                len(data),
                time[0],
                time[-1])               

    
    #Plot data line over time in given color
    return ax.plot(
            time,
            data,
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
                    return 'total population'
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
        
def _safe_fig(title,t_0,t_max,step,family_sizes):
    """Safes the figure with the respective title.
    
    """
    mating_scheme = _get_name_from_mating_scheme()
            
#    initial_state = '({}x{},{})'.format(model['number_of_families'],int(model['K_list'][0]*model['initial_rel_quantity_per_family']),int(model['K_list'][0]*model['initial_rel_mutation_load']))
    path = "{}_{}(K={},n_loci={}".format(title,mating_scheme,model['K_list'],model['number_of_loci'])
    
    if t_0>0:
        path += ',t_0={}'.format(t_0)
#    if t_max != None:
#        path += ',t_max={}'.format(t_max)
    if model['num_runs'] != 0:
        path += ',n_run={}'.format(model['num_runs'])
    if step != 1:
        path += ',step={}'.format(step)
    if family_sizes != None:
        path += ' ,family_sizes={}'.format(family_sizes)
        
    path += ').pdf'
    #safe image 
    plt.savefig(os.path.join(r'/home/larocca/Dokumente/Simulationen/Diploid_Model_Two_Loci/out_analysis/figures',
                                 path))

def _get_name_from_mating_scheme():
    """Returns the mating scheme as a string variable from the global model variable
    
    """
    if model['mating_scheme'][0] == 'dynasty':
        return 'dynasty(nmr={})'.format(model['mating_scheme'][1])
    elif model['mating_scheme'][0] == 'random':
        return 'random'
    elif model['mating_scheme'][0] == 'western':
        return 'western(mr={})'.format(model['mating_scheme'][1])
    elif model['mating_scheme'][0] == 'consang':
        return 'consanguin({},{},{},{})'.format(model['mating_scheme'][1],model['mating_scheme'][2],model['mating_scheme'][3],model['mating_scheme'][4])
    elif model['mating_scheme'][0] == '%in_out_family':
        return '{}%infamily_{}%outfamily'.format(model['mating_scheme'][1],model['mating_scheme'][2])
    else:
        raise ValueError('Mating scheme *{}* in _safe_fig unknown.'.format(model['mating_scheme'][0]))
        
if __name__ == "__main__":
    #Remember to manually set the model name also in plotting_functions.py
    model_name = 'Families'
    
    #Absolute dir the script is in. Needs to be changed manually if the file is moveed
    script_dir = r'/home/larocca/Dokumente/Simulationen/Diploid_Model_Two_Loci'
    
    #load model specific variables such as carrying capacity, initial trait and marker, etc.
    rel_path = "model_specifications/{}.json".format(model_name)
    abs_file_path = os.path.join(script_dir, rel_path)
    model = json.load(open(abs_file_path), encoding="utf-8")
    
    #setup the gene distribution if necessary
    if type(model['gene_distribution']) == list:
        model['number_of_loci'] = model['gene_distribution'][-1]

    parameters_K = [
            (500,10000),
            (1000,15000),
            (2000,20000),
            (5000,25000),
            (10000,30000),
            (20000,35000),
            (50000,70000),
            (100000,100000)]
    
    parameters_n_loci = [50,300,500,1000,1500,2000,2500,3000,5000]
    
#    #set model parameters manually
    model['K_list'] = [500,10000]
    model['t_max'] = 50000
#    model['initial_rel_mutation_load'] = 12
#    model['initial_rel_quantity_per_family'] = 0.0002
#    model['inv_mutation_rate'] = 'per bp'
    
#    model['mating_scheme'] = ['western',0.01]
#    model['mating_scheme'] = ['dynasty',0.1]
#    model['mating_scheme'] = ['random']
    model['mating_scheme'] = ['consang',1,0,0,0]
    
    model['num_runs'] = 5
    model['number_of_loci'] = 1000
#    model['number_of_families_list'] = [1,20]   #[[5,100],[10,200],[15,300],[20,400],[25,500],[50,1000]]
#    family_sizes = [int(K/num_fam) for K,num_fam in zip(model['K_list'],model['number_of_families_list'])]
    family_sizes = [100,100]
    
#----vvvvvv--------------    

    #set the relativ path
#    rel_path = "out_analysis/data/statistic_{}(K={},nloci={},mating={},n_runs={}".format(
#            model_name,
#            model['K_list'],
#            model['number_of_loci'],
#            model['mating_scheme'],
#            model['num_runs']
#            )
#    if family_sizes != None:
#        rel_path += ',family_sizes={}).pickle'.format(family_sizes)
#    else:
#        rel_path += ').pickle'
#
#    # Load locations from pickle file
#    abs_file_path = os.path.join(script_dir, rel_path)
#    with open(abs_file_path, "rb") as in_file:
#        stats_by_time = pickle.load(in_file)

#---Fraction-vvv--------------------

#    data = {}
#    
#    for mating_scheme in [['random'],['consang',1,0,0,0]]:
#        #set the relativ path
#        rel_path = "out_analysis/data/statistic_{}(K={},nloci={},mating={},n_runs={}".format(
#                model_name,
#                model['K_list'],
#                model['number_of_loci'],
#                mating_scheme,
#                model['num_runs']
#                )
#        if family_sizes != None:
#            rel_path += ',family_sizes={}).pickle'.format(family_sizes)
#        else:
#            rel_path += ').pickle'
#            
#        # Load locations from pickle file
#        abs_file_path = os.path.join(script_dir, rel_path)
#        with open(abs_file_path, "rb") as in_file:
#            data[mating_scheme[0]] = pickle.load(in_file)
#            
#    stats_by_time = []
#
#    for t in range(len(data['random'])):
#        ill_random = data['random'][t]['number_false']
#        pop_size_random = sum(data['random'][t]['pop_state'].values())
#        ill_consang = data['consang'][t]['number_false']
#        pop_size_consang = sum(data['consang'][t]['pop_state'].values())
#        
#        if ill_random == 0:
#            nf = 10
#        else:
#            nf = (ill_consang/pop_size_consang)/(ill_random/pop_size_random)
#                        
#        stats_by_time.append({
#                'pop_state': data['consang'][t]['pop_state'],
#                'mutation_load': data['random'][t]['mutation_load'],
##                'number_false': (ill_consang/pop_size_consang) - (ill_random/pop_size_random)
#                'number_false': nf
#
#                })

#----------------------------------------    
                
    #import model specific functions such as birthrate, mutation operator and probabilities, etc.
    spc = importlib.import_module('model_specifications.{}'.format(model_name))
    
    plot_dict_abs = {('total','total'):_plot_abs_quantity}
    
    plot_dict_rel = {
#            ('total',4,'total','total'):_plot_rel_quantity,            
#            ('total',3,'total','total'):_plot_rel_quantity,
#            ('total',2,'total','total'):_plot_rel_quantity,
#            ('total',1,'total','total'):_plot_rel_quantity,
            ('total','>0','total','total'):_plot_rel_quantity,
#            ('homo','total','total','total'):_plot_rel_quantity,
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
#            ('number_false', 'consanguin - random','#fbb900',None, False):_plot_data,
            ('number_false', 'consanguin / random','#2a7e1d',None, False):_plot_data,
#            ('number_of_types','Number of different family-types','g',None): _plot_data,
#            ('mutation_counter', 'Total umber of mutation','g'):_plot_data,
#             ('active_families','# families alive','lightgrey',None):_plot_data,
            }
    
    # key : data_type, label, color, marker, average_line // value : function
    plot_dict_data_rel2 = {
            ('mutation_load','Relative mutation load','#FF0000',None,False):_plot_data_over_pop_size,
#            ('bla','Relative Entropy of family sizes','b',None):_plot_rel_entropy,
            }
    
    plot_dict_data_rel1 = {
            ('number_false','% of diseased individual','#FFA500',None,False):_plot_data_over_pop_size,
#            ('bla','Relative Entropy of family sizes','b',None):_plot_rel_entropy,
            }
#choose title
#    title = 'rel_entropy'
#    title='%_indv_mut_load'
    title='ill_individual_mutation_load_with_rec'
#    title = 'family_types_entropy'
#    title = 'hom_het'
    
    tmax, stepp = None, 1
#    tmax, stepp = 2500, 1
    
#    plot_something_over_time_three_axis(stats_by_time,plot_dict_abs,plot_dict_data_rel1,plot_dict_data_rel2,family_sizes=family_sizes)
#    plot_something_over_time_two_axis(stats_by_time,plot_dict_data,plot_dict_data_rel1,safe=title)
#    plot_something_over_time_one_axis(stats_by_time,plot_dict_data_rel)
#    plot_something_over_time_one_axis(stats_by_time,plot_dict_abs)   
   
#-----------------------------------------------------------------------------------------------------------------------
    
    #or "number_false" or "mutation_load"
#    model['mating_scheme'] = ['random']
#    model['mating_scheme'] = ['consang',1,0,0,0]
#    data_type = 'number_false'
#    averages = _get_averages_different_parameters(parameters_K,parameters_n_loci,data_type)
#    
#    # Store list with locations after each round in a pickle file
#    rel_path = "out_analysis/data/averages_per_pop_size_nloci({}_{})".format(data_type,model['mating_scheme'][0])
#    abs_file_path = os.path.join(script_dir, rel_path)
#    with open(abs_file_path, "wb") as out_file:
#        pickle.dump(averages, out_file)
#
#        
#    #Or "n_loci" or "pop_size"
#    plot_averages(averages,"pop_size",data_type,safe=True)
#    plot_averages(averages,"n_loci",data_type,safe=True)
#    
#    print('Done plotting')
    
    
    
      # Load locations from pickle file
    data_type = "mutation_load"
    model['mating_scheme'][0] = "random"
    rel_path = "out_analysis/data/averages_per_pop_size_nloci({}_{}).pickle".format(data_type,model['mating_scheme'][0])
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "rb") as in_file:
        averages = pickle.load(in_file)
        
    print(averages)
  
    