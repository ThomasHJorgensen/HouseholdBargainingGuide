import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import FormatStrFormatter
from matplotlib.transforms import Bbox

from itertools import product

color_palettes = {
            'IBM': { # Color blind friendly
                'red'       : '#dc267f',
                'blue'      : '#648fff',
                'purple'    : '#785ef0',
                'orange'    : '#fe6100',
                'yellow'    : '#ffbe0b',
                'black'     : '#000000',
            },
            'standard': { # matplotlib standard
                'red'      : '#d62728',
                'blue'     : '#1f77b4',
                'orange'   : '#ff7f0e',
                'green'    : '#2ca02c',
                'purple'   : '#9467bd',
                'brown'    : '#8c564b',
                'pink'     : '#e377c2',
                'gray'     : '#7f7f7f',
                'olive'    : '#bcbd22',
                'cyan'     : '#17becf',
            }
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
        self.output_path = '../figures/'
        self.font_size = 12
        self.subplot_legends = 'all' # all, shared, 1, 2, ..., 16
        self.shared_legend = False
        self.save_as = ''
        self.size = 'A4'
        
        self.color_palettes = {
            'IBM': { # Color blind friendly
                'red'       : '#dc267f',
                'blue'      : '#648fff',
                'purple'    : '#785ef0',
                'orange'    : '#fe6100',
                'yellow'    : '#ffbe0b',
                'black'     : '#000000',
            },
            'standard': { # matplotlib standard
                'red'       : '#d62728',
                'blue'      : '#1f77b4',
                'orange'    : '#ff7f0e',
                'green'     : '#2ca02c',
                'purple'    : '#9467bd',
                'brown'     : '#8c564b',
                'pink'      : '#e377c2',
                'gray'      : '#7f7f7f',
                'olive'     : '#bcbd22',
                'cyan'      : '#17becf',
            },
            'CEBI': { # CEBI standard
                'red'       : '#A7002A',
                'blue'      : '#278094',
                'yellow'    : '#FBBA1D',
                'green'     : '#56882D',
                'dark_blue' : '#0C2C4D',
            },
        }
        self.color_palette = self.color_palettes['CEBI']
        # KU red a31d20, 901a1e
        
        
        
        # self.color_palette = {
        #     'red'           : '#ae2012',
        #     'turquoise'     :'#0a9396',
        #     'light_orange'  : '#ee9b00',
        #     'light_turquoise': '#94d2bd',
        #     'dark_turquise'  : '#005f73',
        #     'orange'        : '#ca6702',
        #     'dark_blue'     : '#001219',
        #     'dark_red'      : '#9b2226',
        #     'dark_orange'   : '#bb3e03',
        #     'baish'         : '#e9d8a6',
        # }
        
        self.linestyle_palette = {'solid'   : '-',
                                  'dashed'  : '--',
                                  'dotted'  : ':',
                                  'dashdot' : '-.',
                                  'longdash': (0, (10,3)),
                                  }
        
        self.marker_palette = {'circle'     : 'o',
                               'square'     : 's',
                               'diamond'    : 'D',
                               'star'       : '*',
                               'pentagon'   : 'P',
                               }
        
        # Apply plot settings
        self.labels = self.infer_labels()
        self.set_curve_settings(colors = self.color_palette, 
                                linestyles= '-', 
                                markers= '')
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
        self.set_x_grid()
        self.set_y_grid()
        if 'subtitles' in kwargs.keys(): 
            self.set_subtitles(kwargs['subtitles'])
        if 'labels' in kwargs.keys():
            self.set_labels(kwargs['labels'])
        if 'x_grid' in kwargs.keys():
            self.set_x_grid(kwargs['x_grid'])
        if 'y_grid' in kwargs.keys():
            self.set_y_grid(kwargs['y_grid'])
            
        self.set_curve_settings(**kwargs)
            
        plt.style.use('default')  # Use default style (no grid lines)
        plt.rcParams.update({'font.size': self.font_size, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'legend.frameon': False, 'lines.linewidth':2.0})  # Set font size and family for all text in figure
    
    def circle_through(self, iterable, length):
        return [iterable[i % len(iterable)] for i in range(length)]
    
    def set_x_grid(self, x_grid=None):
        if x_grid is None:
            self.x_grid = [None] * self.num_subplots
        elif type(x_grid) == list:
            self.x_grid = x_grid
        else:
            self.x_grid = [x_grid] * self.num_subplots
            
    def set_y_grid(self, y_grid=None):
        if y_grid is None:
            self.y_grid = [None] * self.num_subplots
        elif type(y_grid) == list:
            self.y_grid = y_grid
        else:
            self.y_grid = [y_grid] * self.num_subplots
      
    def set_curve_settings(self, **kwargs):
        num_curves = len(self.labels)
        
        for attr in ['colors', 'linestyles', 'markers']:
            if attr in kwargs.keys():
                if type(kwargs[attr]) == str:
                    kwargs[attr] = [kwargs[attr]]
                if type(kwargs[attr]) == dict:
                    kwargs[attr] = list(kwargs[attr].values())
                values = self.circle_through(kwargs[attr], num_curves)
                setattr(self,attr,values)
    
    def plot_color_palette(self):
        # make pie chart of colors
        fig, ax = plt.subplots()
        ax.pie([1]*len(self.color_palette), labels=self.color_palette.keys(), colors=self.color_palette.values())
        
        # close 
        plt.close()
        
        return fig
    
    
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
    
    
    def set_size(self, fig, size='A4'):
        rows, columns = self.dimensions
        if size == 'A4':
            subplot_width = (8.27-2*0.75) / columns
        elif size == 'A3':
            subplot_width = (11.69-2*0.75) / columns
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
        
        # ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        # ax.yaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
        # ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        # ax.yaxis.set_minor_formatter(StrMethodFormatter('{x:,.0f}'))
        # ax.tick_params(axis='both', which='major', labelsize=self.font_size)
        # ax.tick_params(axis='both', which='minor', labelsize=self.font_size)
        
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
        self.set_size(fig, size=self.size)

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
        fig.savefig(self.output_path + figname + self.filetype, dpi=300)
        self.save_subplots(fig, figname)
        
    def save_subplots(self, fig, figname):
        # get axes from fig
        axes = fig.get_axes()
        
        ### Save individual subplots
        for i, ax in enumerate(axes):
            subtitle = self.subtitles[i]
            if subtitle == '':
                subtitle = str(i)
            subtitle = subtitle.replace('\'', '').replace(' ', '_')
            fig.savefig(self.output_path + figname + '_' + subtitle + self.filetype, dpi=300, bbox_inches=self.full_extent(ax, 0.04).transformed(fig.dpi_scale_trans.inverted()))
    
    def add_shared_legend(self, fig, axes):
        lines, labels = axes[0].get_legend_handles_labels()
        fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.05), bbox_transform=fig.transFigure, ncol=3)
    
    def remove_subplot_legends(self, axes):
        for ax in axes:
            ax.get_legend().remove()

        
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
                    x_grid = self.x_grid[i_plot]
                    y_grid = self.y_grid[i_plot]
                    
                    ax = plot_function(model, var, index, ax = ax, title = subtitle, label = label, x_grid = x_grid, y_grid=y_grid, color=self.colors[i_label], linestyle=self.linestyles[i_label], marker=self.markers[i], **subplot_settings)
                    
        if self.subplot_legends == 'shared':
            self.add_shared_legend(fig, axes)
            ax_mask = [ax for i, ax in enumerate(axes.flatten()) if i < self.num_subplots]
            self.remove_subplot_legends(ax_mask)
        elif type(self.subplot_legends) == int:
            ax_mask = [ax for i, ax in enumerate(axes.flatten()) if (i != self.subplot_legends) and i < self.num_subplots]
            self.remove_subplot_legends(ax_mask)
            
                    
        plt.tight_layout()
        
        if self.save_as != '':
            self.save_plot(fig, self.save_as)
        
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



def plot_surplus(model, vars, idx, label='', add_lines=True, title=None, ax=None, **kwargs):
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
    
    
    ax.plot(par.grid_power,Sw, label='Woman', color=color_palettes['IBM']['red'])
    ax.plot(par.grid_power,Sm, label='Man', color=color_palettes['IBM']['blue'])
    ax.plot(par.grid_power,np.zeros(par.num_power), color='black',linestyle='--') # add zero line
    
    #add vertical line at start_power
    if add_lines:
        start_power = par.grid_power[iP]
        ax.axvline(x=start_power, color=color_palettes['standard']['gray'],linestyle='--', label='$\mu_{t-1}$')

        #add vertical line at end_power if it exists
        end_power = power[idx_couple(iP)]
        if end_power >= 0:
            ax.axvline(x=end_power, color=color_palettes['standard']['gray'],linestyle=':', label='$\mu_{t}$')

            # make a one directional arrow from start_power to end_power lower than the x axis
            ax.annotate("", xy=(end_power, -0.1), xytext=(start_power, -0.1), arrowprops=dict(arrowstyle="->"))
    
    # Layout
    ax.set_xlabel('Power of woman')
    ax.set_ylabel('Marriage surplus')
    if title is not None:
        ax.set_title(title)
    ax.legend()
        
    return ax


    
def infer_index(model, variable, index, z):
    
    index_dimension = {
        't': model.par.T,
        'iP': model.par.num_power,
        'iL': model.par.num_love,
        'iA': model.par.num_A
    }
    
    inv_index = {v: k for k, v in index_dimension.items()}
    
    var = getattr
    index_names = [inv_index[i] for i in variable.shape]
    
    z_to_index = {'time': 't', 
                  'power': 'iP', 
                  'love': 'iL', 
                  'assets': 'iA',
                  'assets_woman': 'iA',
                  'assets_man': 'iA',
                  }
    
    idx_slice = [np.nan for i in index_names]
    for i, i_name in enumerate(index_names):
        if i_name == z_to_index[z]:
            idx_slice[i] = slice(None, None, None)
        elif i_name in index_names:
            idx_slice[i] = index[i]
        else:
            print('Warning: index not found')
    
    return tuple(idx_slice)

def get_z_grid(model, z):
    if z == 'time':
        return np.arange(model.par.T)
    if z == 'power':
        return model.par.grid_power
    if z == 'love':
        return model.par.grid_love
    if z == 'assets':
        return model.par.grid_A
    if z == 'assets_woman':
        return model.par.grid_Aw
    if z == 'assets_man':
        return model.par.grid_Am
    
        
def plot_var_over_z(model, variable, index, z, namespace='sol', x_grid=None, y_grid=None, label='', title=None, ax=None, **kwargs):

    nmspc = getattr(model, namespace)
    
    z_grid = get_z_grid(model, z)
    
    if ax is None:
        fig, ax = plt.subplots()
        
    var = getattr(nmspc, variable)
    idx_slice = infer_index(model, var, index, z)
    y = var[idx_slice]
    
    ax.plot(z_grid, y, label=label, **kwargs)
    
    # Layout
    ax.set_xlabel(z.capitalize())
    ax.set_ylabel('')
    if title is not None:
        ax.set_title(title)
    if x_grid is not None:
        ax.set_xlim([x_grid[0], x_grid[-1]])
    if y_grid is not None:
        ax.set_ylim([y_grid[0], y_grid[-1]])
    ax.legend()
    
    return ax

def plot_var_over_assets(*args, **kwargs):
    return plot_var_over_z(*args, z='assets', **kwargs)

def plot_var_over_time(*args, **kwargs):
    return plot_var_over_z(*args, z='time', **kwargs)

def plot_var_over_power(*args, **kwargs):
    return plot_var_over_z(*args, z='power', **kwargs)

def plot_var_over_love(*args, **kwargs):
    return plot_var_over_z(*args, z='love', **kwargs)


def plot_simulated(model, variable, function, label='', title=None, x_grid=None, y_grid=None, ax=None, **kwargs):
    
    # it is possible to add a where condition, but I don't know how to apply the where condition specifically to the each model
    if ax is None:
        fig, ax = plt.subplots()
        
    y = function(getattr(model.sim, variable))
    ax.plot(y, label=label, **kwargs)
    
    # Layout
    ax.set_xlabel('Time')
    ax.set_ylabel('')
    if title is not None:
        ax.set_title(title)
    if x_grid is not None:
        ax.set_xlim([x_grid[0], x_grid[-1]])
    if y_grid is not None:
        ax.set_ylim([y_grid[0], y_grid[-1]])
    ax.legend()
    
    return ax
    

        
