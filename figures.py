import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import FormatStrFormatter
from matplotlib.transforms import Bbox

from itertools import product

import warnings
warning_text_to_ignore = "No artists with labels found to put in legend.  Note that artists whose label start with an underscore are ignored when legend() is called with no argument."
warning_text_to_ignore = "No artists with labels found to put in legend."
warnings.filterwarnings("ignore", message=warning_text_to_ignore)


linestyles = ['-','--','-.',':',':']
markers = ['o','s','D','*','P']
linewidth = 1.5
font_size = 12
font = {'size':font_size}
matplotlib.rc('font', **font)
plt.rcParams.update({'figure.max_open_warning': 0,'text.usetex': False})
path = 'output/'
filetype = '.pdf'

plt.style.use('default')  # Use default style (no grid lines)
plt.rcParams.update({'font.size': font_size, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})  # Set font size and family for all text in figure
    


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

formats = {
    -2: {'font_size': 24},
    1: {'font_size': 16},
    2: {'font_size': 12},
    3: {'font_size': 8},
}

class ModelPlot():
    
    def __init__(self, plot_type, models=None, variables=None, indexes=None, **kwargs):
        
        # set plot type
        self.plot_types = ['models', 'variables', 'indexes']
        self.set_plot_type(plot_type)
        
        # Initialize models
        self.set_models(models)
        
        # Initialize variables
        self.set_variables(variables)
        
        # Initialize indexes
        self.set_indexes(indexes)
        
        # Other
        self.filetype = '.pdf'
        self.path = 'output/'
        self.font_size = 12
        self.shared_legend = False
        self.save = False
        
        # Apply plot settings
        self.set_plot_settings(**kwargs)
        
            
    def set_plot_settings(self, **kwargs):
        
        # Set attributes
        for key, value in kwargs.items():
            if key == 'plot_type':
                self.set_plot_type(value)
            elif key == 'models':
                self.set_models(value)
            elif key == 'variables':
                self.set_variables(value)
            elif key == 'indexes':
                self.set_indexes(value)
            elif key in self.__dict__.keys():
                setattr(self, key, value)
            elif key == 'subplot_settings':
                pass
            else: 
                print(f'Warning: {key} is not a valid attribute')
        
        self.infer_plot()
        if 'subtitles' in kwargs.keys(): 
            self.set_subtitles(kwargs['subtitles'])
        if 'labels' in kwargs.keys():
            self.set_labels(kwargs['labels'])
            
        plt.style.use('default')  # Use default style (no grid lines)
        plt.rcParams.update({'font.size': self.font_size, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})  # Set font size and family for all text in figure
        
    def set_plot_type(self, plot_type):
        if plot_type in self.plot_types:
            self.plot_type = plot_type
        else:
            print('Warning: plot_type should be "models", "variables", or "indexes"')
        
    def set_models(self, models=None):
        self.models = {}
        if models is not None:
            self.add_models(models)
        
    def add_models(self, models):
        if type(models) == dict:
            self.models.update(models)
        elif type(models) == list:
            self.models.update({model.name: model for model in models})
        else:
            try:
                self.models[models.name] = models
            except:
                print('Warning: models should dictionary, list, or model (EconModel)')
            
    def add_variables(self, variables):
        if type(variables) == dict:
            self.variables.update(variables)
        if type(variables) == list:
            self.variables.update({var: var for var in variables})
        elif type(variables) == str:
            self.variables.update({variables: variables})
        else:
            print('Warning: variables should be a dictionary, list, or string')
            
    def set_variables(self, variables=None):
        self.variables = {}
        if variables is not None:
            self.add_variables(variables)
            
    def add_indexes(self, indexes):
        if type(indexes) == dict:
            self.indexes.update(indexes)
        elif type(indexes) == list:
            self.indexes.update({str(idx): idx for idx in indexes})
        elif type(indexes) == tuple:
            self.indexes.update({str(indexes): indexes})
        else:
            print('Warning: indexes should be a dictionary, list, or tuple')
        
    def set_indexes(self, indexes=None):
        self.indexes = {}
        if indexes is not None:
            self.add_indexes(indexes)
    
    def infer_plot(self):
        self.num_subplots = self.infer_number_of_subplots()
        self.dimensions = self.infer_dimensions()
        self.subtitles = self.infer_subtitles()
        self.labels = self.infer_labels()
    
    def infer_number_of_subplots(self):
        return len(getattr(self, self.plot_type))
    
    def infer_dimensions(self):
        num_subplots = self.infer_number_of_subplots()
        if num_subplots == 1:
            dimensions = (1, 1) # one column
        elif num_subplots <= 8:
            dimensions = ((num_subplots+1)//2, 2) # two columns
        elif num_subplots <= 12:
            dimensions = ((num_subplots+2)//3, 3) # three columns
        elif num_subplots <= 16:
            dimensions = ((num_subplots+3)//4, 4) # four columns
        else: 
            print('Warning: Cannot handle more than 16 plots')
            dimensions = (1,1)
                  
        return dimensions
    
    def infer_subtitles(self):
        return [name for name in getattr(self, self.plot_type).keys()] 
    
    def set_subtitles(self, subtitles=None):
        num_subplots = self.infer_number_of_subplots()
        
        if subtitles is None:
            return self.infer_subtitles()
        elif len(subtitles) != num_subplots:
            print('Warning: number of subtitles does not match number of models. Subtitles set to None.')
            return [None] * num_subplots
        else:
            return subtitles
        
    def infer_labels(self):
        label_types = []
        names = [name for name in self.plot_types if name != self.plot_type]
        for i, label_type in enumerate(names):
            labels = [name for name in getattr(self, label_type).keys()]
            if len(labels) > 1:
                label_types.append(labels)
        
        # Use itertools.product to generate the Cartesian product
        cartesian_product = list(product(*label_types))

        # Combine the strings in each tuple
        labels =  [', '.join(combination) for combination in cartesian_product]
        
        return labels
    
    def set_labels(self, labels=None):
        inferred_labels = self.infer_labels()
        num_labels = len(inferred_labels)
        
        if labels is None:
            return self.infer_labels()
        elif len(labels) != num_labels:
            print('Warning: number of subtitles does not match number of models. Subtitles set to None.')
            return [None] * num_labels
        else:
            return labels
    
    
    def set_size(self, fig):
        rows, columns = self.dimensions

        subplot_width = (8.27-2*0.75) / columns
        subplot_height = subplot_width / 1.414
        fig_width = subplot_width * columns
        fig_height = subplot_height * rows
        
        fig.set_size_inches(fig_width, fig_height)
        
    def format_ax(self, ax):
        # Apply font size to tick labels on both axes
        ax.tick_params(axis='both', labelsize=self.font_size)
        
        # Apply font size to axes labels
        ax.xaxis.label.set_fontsize(self.font_size)
        ax.yaxis.label.set_fontsize(self.font_size)
        
        # Apply slightly larger font size to titles
        ax.title.set_fontsize(self.font_size + 2)
                                     
        # Change the font size of y-axis tick labels (including scientific notation)
        ax.yaxis.get_offset_text().set_fontsize(self.font_size - 2)
        
    def create_fig_and_axes(self):
    
        # 1. Setup
        rows, columns = self.dimensions

        # 2. Make figure
        ## a. Make figure and axes
        fig, axes = plt.subplots(rows, columns, squeeze=False)
        axes = axes.flatten()

        ## b. Remove subplots that are not used
        for i in range(self.num_subplots, rows * columns):
            axes[i].axis('off')

        ## c. Apply format to all subplots
        for ax in axes:
            ax = self.format_ax(ax)
            
        # set size
        self.set_size(fig)

        return fig, axes
    
    
        
    def get_plot_iter(self, model_iter, var_iter, index_iter):
        if self.plot_type == 'models':
            return model_iter
        elif self.plot_type == 'variables':
            return var_iter
        elif self.plot_type == 'indexes':
            return index_iter
        
    def get_label_iter(self, model_iter, var_iter, index_iter):
        if self.plot_type == 'models':
            return var_iter + index_iter
        elif self.plot_type == 'variables':
            return model_iter + index_iter
        elif self.plot_type == 'indexes':
            return model_iter + var_iter
        
    def save_plot(self, fig, figname=None):
        fig.savefig(self.path + figname + self.filetype, dpi=300)
        self.save_subplots(fig, figname)
        
    def save_subplots(self, fig, figname):
        # get axes from fig
        axes = fig.get_axes()
        
        ### Save individual subplots
        for i, ax in enumerate(axes):
            subtitle = self.subtitles[i]
            if subtitle == '':
                subtitle = str(i)
            fig.savefig(self.path + figname + '_' + subtitle + self.filetype, dpi=300, bbox_inches=full_extent(ax, 0.04).transform)
    
    def add_shared_legend(self, fig, axes):
        lines, labels = axes[0].get_legend_handles_labels()
        fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.05), bbox_transform=fig.transFigure, ncol=3)
    
    def remove_subplot_legends(self, axes):
        for i in range(self.num_subplots):
            axes[i].get_legend().remove()

        
    def plot_using(self, plot_function, **kwargs):
        
        self.set_plot_settings(**kwargs)
        if 'subplot_settings' in kwargs:
            subplot_settings = kwargs['subplot_settings']
        else:
            subplot_settings = {}
        
        fig, axes = self.create_fig_and_axes()

        for m, model in enumerate(self.models.values()):
            for v, var in enumerate(self.variables.values()):
                for i, index in enumerate(self.indexes.values()):
                    i_plot = self.get_plot_iter(m, v, i)
                    i_label = self.get_label_iter(m, v, i)
                    subtitle = self.subtitles[i_plot]
                    label = self.labels[i_label]
                    ax = axes[i_plot]
                    
                    ax = plot_function(model, var, index, ax = ax, title = subtitle, label = label, **subplot_settings)
                    
        if self.shared_legend:
            self.add_shared_legend(fig, axes)
            self.remove_subplot_legends(axes)
                    
        plt.tight_layout()
        
        if self.save:
            self.save_plot(fig, self.figname)
        
        plt.close()
        
        return fig


    
    
    def full_extent(self, ax, pad=0.0):
        """Get the full extent of an axes, including axes labels, tick labels, and
        titles."""
        # For text objects, we need to draw the figure first, otherwise the extents
        # are undefined.
        ax.figure.canvas.draw()
        items = []
        # items += ax.get_xticklabels() + ax.get_yticklabels() 
        items += [ax, ax.xaxis.label, ax.yaxis.label]
        # items += [ax, ax.title]
        
        # items += [ax.get_xaxis().get_label(), ax.get_yaxis().get_label()]
        bbox = Bbox.union([item.get_window_extent() for item in items])

        # return bbox.padded(4)
        return bbox.expanded(1.0 + pad, 1.0 + pad)

    
def plot_surplus(model, vars, idx, label='', add_lines=True, title=None, ax=None):
    par = model.par
    sol = model.sol
    power = sol.power
    
    # Create indices
    try:
        t, iP, iL, iA = idx
    except:
        print('Warning: idx not for couple')
        return
    idx_couple = lambda iP: (t,iP,iL,iA)
    idx_single = (t,iA)
    
    # Calculate surplus
    Sw = sol.Vw_couple_to_couple[t,:,iL,iA]-sol.Vw_couple_to_single[idx_single]
    Sm = sol.Vm_couple_to_couple[t,:,iL,iA]-sol.Vm_couple_to_single[idx_single]
    
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


    
def infer_index(model, variable, index, x):
    
    index_dimension = {
        't': model.par.T,
        'iP': model.par.num_power,
        'iL': model.par.num_love,
        'iA': model.par.num_A
    }
    
    inv_index = {v: k for k, v in index_dimension.items()}
    
    var = getattr
    index_names = [inv_index[i] for i in variable.shape]
    
    x_to_index = {'time': 't', 
                  'power': 'iP', 
                  'love': 'iL', 
                  'assets': 'iA',
                  'assets_woman': 'iA',
                  'assets_man': 'iA',
                  }
    
    idx_slice = [np.nan for i in index_names]
    for i, i_name in enumerate(index_names):
        if i_name == x_to_index[x]:
            idx_slice[i] = slice(None, None, None)
        elif i_name in index_names:
            idx_slice[i] = index[i]
        else:
            print('Warning: index not found')
    
    return tuple(idx_slice)

def get_x_grid(model, x):
    if x == 'time':
        return np.arange(model.par.T)
    if x == 'power':
        return model.par.grid_power
    if x == 'love':
        return model.par.grid_love
    if x == 'assets':
        return model.par.grid_A
    if x == 'assets_woman':
        return model.par.grid_Aw
    if x == 'assets_man':
        return model.par.grid_Am
    
        
def plot_var_over_x(model, variable, index, x, namespace='sol', label='', title=None, ax=None):

    nmspc = getattr(model, namespace)
    
    x_grid = get_x_grid(model, x)
    
    if ax is None:
        fig, ax = plt.subplots()
        
    var = getattr(nmspc, variable)
    idx_slice = infer_index(model, var, index, x)
    y = var[idx_slice]
    
    ax.plot(x_grid, y, label=label, alpha=0.5)
    
    # Layout
    ax.set_xlabel(x.capitalize())
    ax.set_ylabel('')
    if title is not None:
        ax.set_title(title)
    ax.legend()
    
    return ax

def plot_var_over_assets(*args, **kwargs):
    return plot_var_over_x(*args, x='assets', **kwargs)

def plot_var_over_time(*args, **kwargs):
    return plot_var_over_x(*args, x='time', **kwargs)

def plot_var_over_power(*args, **kwargs):
    return plot_var_over_x(*args, x='power', **kwargs)

def plot_var_over_love(*args, **kwargs):
    return plot_var_over_x(*args, x='love', **kwargs)


def plot_simulated(model, variable, function, label='', title=None, ax=None,):
    
    # it is possible to add a where condition, but I don't know how to apply the where condition specifically to the each model
    if ax is None:
        fig, ax = plt.subplots()
        
    y = function(getattr(model.sim, variable))
    ax.plot(y, label=label)
    
    # Layout
    ax.set_xlabel('Time')
    ax.set_ylabel('')
    if title is not None:
        ax.set_title(title)
    ax.legend()
    
    return ax
    

        
