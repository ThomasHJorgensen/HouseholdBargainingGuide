import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

linestyles = ['-','--','-.',':',':']
markers = ['o','s','D','*','P']
linewidth = 1.5
font_size = 18
font = {'size':font_size}
matplotlib.rc('font', **font)
plt.rcParams.update({'figure.max_open_warning': 0,'text.usetex': False})
path = 'output/'
filetype = '.pdf'


colors = {
    'red'         : '#FF5733',
    'red_dark'    : '#AA0000',
    'green'       : '#33FF57',
    'green_dark'  : '#008800',
    'blue'        : '#3366FF',
    'blue_dark'   : '#0033CC',
    'orange'      : '#FFA533',
    'orange_dark' : '#FF6600',
    'purple'      : '#A733FF',
    'purple_dark' : '#6600CC',
    'pink'        : '#FF33B8',
    'pink_dark'   : '#FF0066',
    'brown'       : '#A36633',
    'brown_dark'  : '#663300',
    'gray'        : '#A6A6A6',
    'gray_dark'   : '#333333',
    'cyan'        : '#33FFFF',
    'cyan_dark'   : '#00CCCC',
    'magenta'     : '#FF33FF',
    'magenta_dark': '#FF00CC',
    'black'       : '#000000',
}

# To do. 
# Format scientific notation in top of y axis to 10^x instead of 1e^x in make_fig
# Make make_fig able to fit to landscape view
# Consider how to handle diffs in boolean arrays (e.g. couple)
# Add function to add multiple models to a single plot


def make_fig(num_plots: int, dimensions: tuple):
    rows, columns = dimensions

    # Calculate subplot size and figure size
    subplot_width = (8.27-2*0.75) / columns  # Width of an A4 paper
    subplot_height = subplot_width / 1.414  # Height adjusted to fit width
    fig_width = subplot_width * columns
    fig_height = subplot_height * rows

    # Create a figure with the specified size
    fig, axes = plt.subplots(rows, columns, figsize=(fig_width, fig_height), squeeze=False)
    axes = axes.flatten()

    # Calculate font size
    # font_size = 20 * (subplot_width / 8.27)  # Adjust font size
    
    # Apply style and font size
    plt.style.use('default')  # Use default style (no grid lines)
    plt.rcParams.update({'font.size': font_size, 'font.family': 'serif', 'font.serif': 'Times New Roman'})

    # Remove ticks and labels from extra subplots
    for i in range(num_plots, rows * columns):
        axes[i].axis('off')
        
    # Remove ticks and labels from all subplots
    for ax in axes:
        # Apply font size to tick labels on both axes
        ax.tick_params(axis='both', labelsize=font_size)
        
        # Apply font size to axes labels
        ax.xaxis.label.set_fontsize(font_size)
        ax.yaxis.label.set_fontsize(font_size)
        
        # Apply slightly larger font size to titles
        ax.title.set_fontsize(font_size + 2)
                                     
        # Change the font size of y-axis tick labels (including scientific notation)
        ax.yaxis.get_offset_text().set_fontsize(font_size - 2)  # Adjust the labelsize as needed
        
        # # Change the format of scientific notation on the y-axis to show only 4 digits
        # ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.4g}"))

    return fig, axes
       
    
def model_plot(models, plot_function, *args, subtitles=None, num_plots=None, dim=None, save=False, figname=None, display=False, **kwargs):
    ### 1. Setup
    # a. Place models in a list if they are not already
    if type(models) is dict:
        models = [*models.values()]
    if type(models) is not list:
        models = [models]
        
    # b. Set default values for num_plots and dim
    if num_plots is None:
        num_plots = len(models)
    if dim is None:
        if num_plots == 1: # If there are less than 3 plots, make a single row
            dim = (1, 1)
        else:
            dim = ((num_plots+1)//2, 2)
        if num_plots > 10:
            print('Warning: more than 10 plots. Specify dimensions manually.')
            
    # c. Handle subtitles
    if subtitles is None:
        subtitles = [None] * len(models)
    elif subtitles == 'model_names':
        subtitles = [model.name for model in models]
    elif len(subtitles) != len(models):
        print('Warning: number of subtitles does not match number of models. Subtitles set to None.')
        subtitles = [None] * len(models)
    else:
        subtitles = subtitles
    
    ### 2. Create figure
    # a. initiate figure
    fig, ax = make_fig(num_plots, dim)
    
    # b. Create subplots
    for i, model in enumerate(models):
        subtitle = subtitles[i]
        plot_function(model, *args, ax = ax[i], title=subtitle, **kwargs)
        
    # c. Set tight layout
    plt.tight_layout() 
    
    # d. Save the figure
    if save:
        if figname == None:
            figname = plot_function.__name__
        plt.savefig(path + figname + filetype, dpi=300)
    
    # e. Display the figure
    if display:
        plt.show()
    plt.close()
    
    return fig
    
    
def plot_surplus(model, t, iP, iL, iA, add_lines=True, title=None, ax=None):
    par = model.par
    sol = model.sol
    power = sol.power
    
    # Create indices
    idx_couple = lambda iP: (t,iP,iL,iA)
    idx_single = (t,iA)
    
    # Calculate surplus
    Sw = sol.Vw_remain_couple[t,:,iL,iA]-sol.Vw_single[idx_single]
    Sm = sol.Vm_remain_couple[t,:,iL,iA]-sol.Vm_single[idx_single]
    
    # If ax is not provided, create a new figure and axis
    if ax is None:
        fig, ax = plt.subplots()
    
    
    ax.plot(par.grid_power,Sw, label='Woman', color=colors['red'])
    ax.plot(par.grid_power,Sm, label='Man', color=colors['blue'])
    ax.plot(par.grid_power,np.zeros(par.num_power), color='black',linestyle='--') # add zero line
    
    #add vertical line at start_power
    if add_lines:
        start_power = par.grid_power[iP]
        ax.axvline(x=start_power, color=colors['gray'],linestyle='--', label='$\mu_{t-1}$')

        #add vertical line at end_power if it exists
        end_power = power[idx_couple(iP)]
        if end_power >= 0:
            ax.axvline(x=end_power, color=colors['gray_dark'],linestyle=':', label='$\mu_{t}$')

            # make a one directional arrow from start_power to end_power lower than the x axis
            ax.annotate("", xy=(end_power, -0.1), xytext=(start_power, -0.1), arrowprops=dict(arrowstyle="->"))
    
    # Layout
    ax.set_xlabel('Power of woman')
    ax.set_ylabel('Marriage surplus')
    if title is not None:
        ax.set_title(title)
    ax.legend()
        
    return ax

def plot_var_over_assets(model, vars, idx, grid_points=False, title=None, ax=None):
    par = model.par
    sol = model.sol
    
    try:
        t, iP, iL, iA = idx
    except:
        t, iA = idx
    
    if grid_points:
        assets = np.arange(par.num_A)
    else:
        assets = par.grid_A

    # If ax is not provided, create a new figure and axis
    if ax is None:
        fig, ax = plt.subplots()
    
    for var in vars:
        y = getattr(sol, var)
        try:
            y = y[t,iP,iL,:]
        except:
            y = y[t,:]
        
        ax.plot(assets, y, label=var, alpha=0.5)
    

    # Layout
    ax.set_xlabel('Assets')
    ax.set_ylabel('')
    if title is not None:
        ax.set_title(title)
    ax.legend()
        
    return ax

def plot_var_over_time(model, vars, idx, title=None, ax=None):
    par = model.par
    sol = model.sol
    
    try:
        t, iP, iL, iA = idx
    except:
        t, iA = idx

    time = np.arange(par.T)

    # If ax is not provided, create a new figure and axis
    if ax is None:
        fig, ax = plt.subplots()
    
    for var in vars:
        y = getattr(sol, var)
        try:
            y = y[:,iP,iL,iA]
        except:
            y = y[:,iA]
        
        ax.plot(time, y, label=var, alpha=0.5)
    

    # Layout
    ax.set_xlabel('Time')
    ax.set_ylabel('')
    if title is not None:
        ax.set_title(title)
    ax.legend()
        
    return ax              
        
def plot_life_cycle(model_list, fig_name):

    # Create a figure with subplots
    fig, ax = make_fig(num_plots=9, dimensions=(3, 3))
    var_list = ('Cw_priv','Cm_priv','Cw_pub','C_tot','A','power','power_idx','love','couple')
    
    for i, var in enumerate(var_list):
        for j, model in enumerate(model_list.values()):
            # Pick out couples (if not the share of couples is plotted)
            if var == 'couple':
                nan = 0.0
            else:
                I = model.sim.couple < 1
                nan = np.zeros(I.shape)
                nan[I] = np.nan

            # Pick the relevant variable for couples
            y = getattr(model.sim, var)        
            y = np.nanmean(y + nan, axis=0)

            ax[i].plot(y, marker=markers[j], linestyle=linestyles[j], linewidth=linewidth, label=model.spec['latexname'])
            ax[i].set_xlabel('age')
            ax[i].set_title(f'{var}')
            ax[i].legend()
    
    # Set tight layout
    plt.tight_layout()
    
    # Save and display the figure
    plt.savefig(path + 'life_cycle' + fig_name + filetype, dpi=300)
    plt.show()
    plt.close()
        
        

        
