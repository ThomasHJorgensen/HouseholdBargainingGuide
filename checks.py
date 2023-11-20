import numpy as np
# import simplenamespace
from types import SimpleNamespace

#make a simple namespace called par
par = SimpleNamespace()
# Make a loop and add one=1, two=2, three=3, four=4, five=5 to par
for i in range(1,6):
    setattr(par, f'var{i}', i)
    

def difference_in_namespace(namespace_1, namespace_2, relative=False, output='', time=None):
    """
    The function `difference_in_namespace` compares two namespaces and returns a new namespace
    containing the differences between the variables in the two namespaces. Note that the two
    namespaces should have the same variables.
    
    Args:
      namespace_1: The first namespace object to compare.
      namespace_2: The `namespace_2` parameter is the namespace object that you want to compare with
    `namespace_1`. It contains variables that you want to compare with the variables in `namespace_1`.
      relative: The `relative` parameter is a boolean flag that determines whether the differences
    between variables should be calculated as absolute differences (`False`) or relative differences
    (`True`). Defaults to False
    
    Returns:
      a namespace object that contains the differences between the variables in `namespace_1` and
    `namespace_2`.
    """

    # Initialize a new namespace for differences
    namespace_diff = SimpleNamespace()
    
    # Find differences in variables
    for name in namespace_1.__dict__.keys():
        if hasattr(namespace_2, name):
            if time == None:
                var1 = getattr(namespace_1, name)
                var2 = getattr(namespace_2, name)
            else:
                var1 = getattr(namespace_1, name)[time]
                var2 = getattr(namespace_2, name)[time]
            
            # Check if var1 and var2 are boolean arrays
            if np.isnan(var1).any() or np.isnan(var2).any():
                print(f"Variable '{name}' contains nan values and is skipped")
                continue
            if np.issubdtype(var1.dtype, np.bool_) and np.issubdtype(var2.dtype, np.bool_):
                diff = var2 ^ var1  # Use ^ operator for boolean arrays
            if relative:
                diff = np.where(var1 == 0, np.nan, (var2 - var1)/var1)
            else:
                diff = var2 - var1  # Use - operator for other types of arrays
            
            if output == '':
                diff = diff
            elif output == 'max_abs_value':
                diff = max_absolute_value(diff)
            elif output == 'max_abs_index':
                diff = index_of_max_absolute_value(diff, name)
                
            setattr(namespace_diff, name, diff)
        else:
            print(f"Variable '{name}' not found in model2.sim")
            setattr(namespace_diff, name, None)

    return namespace_diff


def max_absolute_value(variable: np.ndarray):
    return np.nanmax(np.abs(variable))
        
def index_of_max_absolute_value(variable: np.ndarray, name):
    # get the index of the maximum absolute value. The function should return a tuple of integers
    # that can be used to index the array
    # index_in_arrays = np.nanargmax(np.abs(variable))
    index_in_arrays = np.where(np.abs(variable) == max_absolute_value(variable))
    if len(index_in_arrays[0]) > 1:
        print(f"Multiple indices for variable: {name}. Return nan")
        # index = tuple(i[0][0] for i in index_in_arrays)
        index = np.nan
    elif len(index_in_arrays[0]) == 0:
        print(f"No indices for variable: {name}. Return nan")
        index = np.nan
    else:
        index = tuple(i[0] for i in index_in_arrays)
    return index
    
    # return np.where(np.abs(variable) == max_absolute_value(variable))

def print_namespace(namespace: SimpleNamespace):
    for name in namespace.__dict__.keys():
        print(f"{name:24}: {getattr(namespace, name)}")
        
def print_specs_table(modelspecs):
    # find the longest name in specs
    max_length = 0
    for key in modelspecs.keys():
        if len(key)>max_length:
            max_length = len(key)
    for key in list(modelspecs.values())[0]['par'].keys():
        if len(key)>max_length:
            max_length = len(key)
            
    # get all keys in par dictionaries
    all_keys = set()
    for spec in modelspecs.values():
        all_keys |= set(spec['par'].keys())
            
    # print header with model names
    print(f"{' '*max_length}|",end='')
    for name in modelspecs.keys():
        print(f"{name.center(max_length)}|",end='')
    print('')
    print('-'*(max_length+1)*(len(modelspecs.keys())+1))
    
    # print parameters
    for p in all_keys:
        print(f"{p.ljust(max_length)}|",end='')
        for name in modelspecs.keys():
            if p in modelspecs[name]['par'].keys(): # use value from modelspecs
                par_value = modelspecs[name]['par'][p]
            else: # use default value
                model = model = brg.HouseholdModelClass()
                par_value = getattr(model.par,p)
            
            print(f"{str(par_value).center(max_length)}|",end='')
        print('')