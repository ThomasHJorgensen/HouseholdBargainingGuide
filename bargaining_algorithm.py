import numpy as np
import numba as nb
import scipy.optimize as optimize

from EconModel import EconModelClass
from consav.grids import nonlinspace
from consav import linear_interp, linear_interp_1d
from consav import quadrature

import time

def check_participation_constraints(power_idx,power,Sw,Sm,idx_single,idx_couple,list_start_as_couple,list_remain_couple,list_trans_to_single, par):

    min_w = Sw[0]
    max_w = Sw[-1]
    min_m = Sm[-1]
    max_m = Sm[0]

    cross_w = (min_w < 0.0) & (max_w > 0.0)
    cross_m = (min_m < 0.0) & (max_m > 0.0)

    always_happy_w = (min_w > 0.0)
    always_happy_m = (min_m > 0.0)

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
            divorce(iP, power_idx, power,idx_single,idx_couple, list_start_as_couple,list_trans_to_single, par)

    # 1c: check if husband is always happy, wife has indifference point
    elif cross_w & always_happy_m: # husband is always happy, wife has indifference point
        # find wife's indifference point
        Low_w = find_left_point(Sw, par) + 1 # lowest point where she has positive surplus
        power_at_zero_w = linear_interp_1d._interp_1d(Sw,par.grid_power,0.0,Low_w-1) # interpolated power at indifference point

        # update case 1c
        for iP in range(par.num_power):
            idx = idx_couple(iP)

            if iP == 0:
                update_to_indifference_point(iP, Low_w, Low_w-1, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par)
            elif iP < Low_w: # update to wife's indifference point
                update_to_indifference_point(iP, Low_w, Low_w-1, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=0)
            else:
                remain(iP, power_idx, power,idx_couple, list_start_as_couple,list_remain_couple, par)
            # if iP < Low_w: # update to wife's indifference point
            #     power_idx[idx] = Low_w
            #     power[idx] = power_at_zero_w
            # else:
            #     power_idx[idx] = iP
            #     power[idx] = par.grid_power[iP]

            # for i,key in enumerate(list_start_as_couple):
            #     if iP == 0:
            #         list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i],power_at_zero_w,Low_w-1)
            #     elif iP < Low_w: # update to wife's indifference point
            #         list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple(0)]    
            #     else: # no bargaining
            #         list_start_as_couple[i][idx] = list_remain_couple[i][iP]                


    # 1d: check if wife is always happy, husband has indifference point
    elif cross_m & always_happy_w: # wife is always happy, husband has indifference point
        # find husband's indifference point
        Low_m = find_left_point(Sm, par) # lowest point where he has positive surplus
        power_at_zero_m = linear_interp_1d._interp_1d(Sm,par.grid_power,0.0,Low_m) # interpolated power at indifference point

        # update case 1d
        for iP in range(par.num_power):
            idx = idx_couple(iP)

            if iP <= Low_m:
                remain(iP, power_idx, power,idx_couple, list_start_as_couple,list_remain_couple, par)
            elif iP == Low_m+1: # update to husbands indifference point
                update_to_indifference_point(iP, Low_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par)
            else:
                update_to_indifference_point(iP, Low_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=Low_m+1)

            # idx = idx_couple(iP)

            # if iP > Low_m:
            #     power_idx[idx] = Low_m
            #     power[idx] = power_at_zero_m
            # else:
            #     power_idx[idx] = iP
            #     power[idx] = par.grid_power[iP]
            
            # for i,key in enumerate(list_start_as_couple):
            #     if iP <= Low_m: # no bargaining
            #         list_start_as_couple[i][idx] = list_remain_couple[i][iP] # no bargaining
            #     elif iP == Low_m+1: # update to husbands indifference point
            #         list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i],power_at_zero_m,Low_m)
            #     else:
            #         list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple(Low_m+1)]


    # 1e: Both have indifference points
    else:
        # find indifference points
        Low_w = find_left_point(Sw, par) + 1 # lowest point where she has positive surplus
        power_at_zero_w = linear_interp_1d._interp_1d(Sw,par.grid_power,0.0,Low_w-1) # interpolated power at indifference point

        Low_m = find_left_point(Sm, par) # lowest point where he has positive surplus
        power_at_zero_m = linear_interp_1d._interp_1d(Sm,par.grid_power,0.0,Low_m) # interpolated power at indifference point

        # step 3: update bargaining power and consumption allocations
        if power_at_zero_w>power_at_zero_m: # no room for bargaining
            for iP in range(par.num_power):
                divorce(iP,power_idx, power,idx_single,idx_couple, list_start_as_couple,list_trans_to_single, par)
        
        else: # bargaining
            for iP in range(par.num_power): 
                idx = idx_couple(iP)

                if iP == 0: #update to woman's indifference point
                    update_to_indifference_point(iP, Low_w, Low_w-1, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par)
                elif iP < Low_w: #re-use precomputed values
                    update_to_indifference_point(iP, Low_w, Low_w-1, power_at_zero_w, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=0)
                
                elif (iP >= Low_w) & (iP <= Low_m): # No change between Low_w and Low_m
                    remain(iP, power_idx, power,idx_couple, list_start_as_couple,list_remain_couple, par)
                
                elif iP == Low_m+1: #update to man's indifference point
                    update_to_indifference_point(iP, Low_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par)
                else: # re-use precomputed values
                    update_to_indifference_point(iP, Low_m, Low_m, power_at_zero_m, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=Low_m+1)
               
                # # update bargaining power
                # if iP<Low_w: # update to woman's indifference point
                #     power_idx[idx] = Low_w
                #     power[idx] = power_at_zero_w
                # elif iP>Low_m: # update to man's indifference point
                #     power_idx[idx] = Low_m
                #     power[idx] = power_at_zero_m
                # else: # bargaining power constant between Low_w and Low_m
                #     power_idx[idx] = iP
                #     power[idx] = par.grid_power[iP]

                # # update consumption allocations
                # for i,key in enumerate(list_start_as_couple):
                #     if (iP==0): # update to woman's indifference point
                #         list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i],power_at_zero_w,Low_w-1)
                #     elif iP < Low_w: # bargaining power constant until Low_w
                #         list_start_as_couple[i][idx]=list_start_as_couple[i][idx_couple(0)]
                #     elif (iP >= Low_w) & (iP <= Low_m): # No change between Low_w and Low_m   
                #         list_start_as_couple[i][idx] = list_remain_couple[i][iP]
                #     elif iP == Low_m+1: # update to man's indifference point
                #         list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i],power_at_zero_m,Low_m)
                #     else:   # bargaining power constant after Low_m
                #         list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple(Low_m+1)]; # re-use that the interpolated values are identical





def divorce(iP,power_idx, power,idx_single,idx_couple, list_start_as_couple,list_trans_to_single, par):
    # overwrite output for couple
    idx = idx_couple(iP)
    for i,key in enumerate(list_start_as_couple):
        list_start_as_couple[i][idx] = list_trans_to_single[i][idx_single]

    power_idx[idx] = -1
    power[idx] = -1.0


def remain(iP,power_idx, power,idx_couple, list_start_as_couple,list_remain_couple, par):
    idx = idx_couple(iP)
    power_idx[idx] = iP
    power[idx] = par.grid_power[iP]

    for i,key in enumerate(list_start_as_couple):
        list_start_as_couple[i][idx] = list_remain_couple[i][iP]


def update_to_indifference_point(iP, Low, left_point, power_at_zero, power_idx, power, idx_couple, list_start_as_couple, list_remain_couple, par, sol_idx=None):
    idx = idx_couple(iP)
    power_idx[idx] = Low
    power[idx] = power_at_zero

    # update solution arrays
    if sol_idx is None: # pre-computation not done
        for i,key in enumerate(list_start_as_couple):
            list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i],power_at_zero,left_point)
    else: # pre-computation done - get solution at sol_idx
        for i,key in enumerate(list_start_as_couple):
            list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple(sol_idx)]
        

def find_left_point(S, par):
    # flip if descending
    if S[0] > S[-1]:
        S = -S
    return linear_interp.binary_search(0,len(S), S, 0.0)

