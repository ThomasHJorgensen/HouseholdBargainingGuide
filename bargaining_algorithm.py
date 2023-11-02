import numpy as np
import numba as nb
import scipy.optimize as optimize

from EconModel import EconModelClass
from consav.grids import nonlinspace
from consav import linear_interp, linear_interp_1d
from consav import quadrature

import time

def check_participation_constraints(power_idx, power, Sw, Sm, idx_couple, list_start_as_couple, list_remain_couple, list_trans_to_single, par):
    """
    Checks the participation constraints for each couple and updates the couple status, bargaining power, and consumption allocations accordingly.

    Args:
    - power_idx (list): list of power indices for each couple
    - power (list): list of power values for each couple
    - Sw (list): list of wife's surplus values for each power index
    - Sm (list): list of husband's surplus values for each power index
    - idx_couple (list): list of indices for coupled individuals
    - list_start_as_couple (list): list of indices for individuals who start as a couple
    - list_remain_couple (list): list of indices for individuals who remain as a couple
    - list_trans_to_single (list): list of indices for individuals who transition from a couple to single
    - par (object): object containing model parameters

    Returns:
    - None
    """

    # step 0: Identify key indicators for each couple
    # 0a: min and max surplus for each couple
    min_w = Sw[0]
    max_w = Sw[-1]
    min_m = Sm[-1]
    max_m = Sm[0]

    # 0b: check if wife and husband have indifference points
    cross_w = (min_w < 0.0) & (max_w > 0.0)
    cross_m = (min_m < 0.0) & (max_m > 0.0)

    # 0c: check if wife and husband are always happy
    always_happy_w = (min_w > 0.0)
    always_happy_m = (min_m > 0.0)

    # 0d: check if wife and husband are never happy
    never_happy_w = (max_w < 0.0)
    never_happy_m = (max_m < 0.0)

    # step 1: check end points
    # 1a: check if all values are consistent with marriage
    if always_happy_w & always_happy_m: # all values are consistent with marriage
        for iP in range(par.num_power):
            remain(iP, power_idx, power,idx_couple, list_start_as_couple,list_remain_couple, par)

    # 1b: check if all values are consistent with divorce
    elif never_happy_w | never_happy_m: # no value is consistent with marriage
        for iP in range(par.num_power):
            divorce(iP, power_idx, power,idx_couple, list_start_as_couple,list_trans_to_single)

    # 1c: check if husband is always happy, wife has indifference point
    elif cross_w & always_happy_m: # husband is always happy, wife has indifference point
        # find wife's indifference point
        left_w = find_left_point(Sw) # index left of indifference point
        Low_w = left_w + 1 # lowest point where she has positive surplus
        power_at_zero_w = linear_interp_1d._interp_1d(Sw,par.grid_power,0.0,Low_w-1) # interpolated power at indifference point

        # update case 1c
        for iP in range(par.num_power):
            if iP == 0:
                update_to_indifference_point(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par)
            elif iP < Low_w: # update to wife's indifference point
                update_to_indifference_point(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=0)
            else:
                remain(iP, power_idx, power,idx_couple, list_start_as_couple,list_remain_couple, par)         


    # 1d: check if wife is always happy, husband has indifference point
    elif cross_m & always_happy_w: # wife is always happy, husband has indifference point
        # find husband's indifference point
        left_m = find_left_point(Sm) # index left of indifference point
        Low_m = left_m # lowest point where he has positive surplus
        power_at_zero_m = linear_interp_1d._interp_1d(Sm,par.grid_power,0.0,Low_m) # interpolated power at indifference point

        # update case 1d
        for iP in range(par.num_power):
            if iP <= Low_m:
                remain(iP, power_idx, power,idx_couple, list_start_as_couple,list_remain_couple, par)
            elif iP == Low_m+1: # update to husbands indifference point
                update_to_indifference_point(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par)
            else:
                update_to_indifference_point(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=Low_m+1)


    # 1e: Both have indifference points
    else:
        # find indifference points
        left_w = find_left_point(Sw) # index left of indifference point
        Low_w = left_w + 1 # lowest point where she has positive surplus
        power_at_zero_w = linear_interp_1d._interp_1d(Sw,par.grid_power,0.0,Low_w-1) # interpolated power at indifference point

        left_m = find_left_point(Sm) # index left of indifference point
        Low_m = left_m # lowest point where he has positive surplus
        power_at_zero_m = linear_interp_1d._interp_1d(Sm,par.grid_power,0.0,Low_m) # interpolated power at indifference point

        # step 3: update bargaining power and consumption allocations
        if power_at_zero_w>power_at_zero_m: # no room for bargaining
            for iP in range(par.num_power):
                divorce(iP,power_idx, power,idx_couple, list_start_as_couple,list_trans_to_single)
        
        else: # bargaining
            for iP in range(par.num_power): 
                if iP == 0: #update to woman's indifference point
                    update_to_indifference_point(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par)
                    
                elif iP < Low_w: #re-use precomputed values
                    update_to_indifference_point(iP, left_w, Low_w, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=0)
                    
                elif (iP >= Low_w) & (iP <= Low_m): # No change between Low_w and Low_m
                    remain(iP, power_idx, power,idx_couple, list_start_as_couple,list_remain_couple, par)
                    
                elif iP == Low_m+1: #update to man's indifference point
                    update_to_indifference_point(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par)
                
                else: # re-use precomputed values
                    update_to_indifference_point(iP, left_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=Low_m+1)
               

def divorce(iP,power_idx, power,idx_couple, list_start_as_couple,list_trans_to_single):
    idx = idx_couple(iP)
    power_idx[idx] = -1
    power[idx] = -1.0

    for i in range(len(list_start_as_couple)):
        list_start_as_couple[i][idx] = list_trans_to_single[i]

    power_idx[idx] = -1
    power[idx] = -1.0


def remain(iP,power_idx, power,idx_couple, list_start_as_couple,list_remain_couple, par):
    idx = idx_couple(iP)
    power_idx[idx] = iP
    power[idx] = par.grid_power[iP]

    for i in range(len(list_start_as_couple)):
        list_start_as_couple[i][idx] = list_remain_couple[i][idx]



def update_to_indifference_point(iP, left_point, low_point, power_at_zero, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=None):
    idx = idx_couple(iP)
    power_idx[idx] = low_point
    power[idx] = power_at_zero
    idxgrid_power_couple = idx_couple(np.arange(par.num_power))

    # update solution arrays
    if sol_idx is None: # pre-computation not done
        for i,key in enumerate(list_start_as_couple):
            list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i][idxgrid_power_couple],power_at_zero,left_point)
    else: # pre-computation done - get solution at sol_idx
        for i,key in enumerate(list_start_as_couple):
            list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple(sol_idx)]
        

def find_left_point(S):
    """
    Finds the left point of a given sequence S using binary search.

    Args:
    S (list): A list of numbers representing a sequence.

    Returns:
    int: The index of the left point of the sequence S.
    """
    if S[0] <= S[-1]:
        return linear_interp.binary_search(0,len(S), S, 0.0)
    else:
        return binary_search_over_descending_function(0,len(S), S, 0.0)

def binary_search_over_descending_function(idx,Nx,x,target):
        
    # a. checks
    if target >= x[0]:
        return 0
    elif target <= x[Nx-2]:
        return Nx-2
    
    # b. binary search
    half = Nx//2
    while half:
        imid = idx + half
        if x[imid] >= target:
            idx = imid
        Nx -= half
        half = Nx//2
        
    return idx