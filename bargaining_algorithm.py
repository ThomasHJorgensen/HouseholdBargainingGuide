import numpy as np
import numba as nb
import scipy.optimize as optimize

from EconModel import EconModelClass
from consav.grids import nonlinspace
from consav import linear_interp, linear_interp_1d
from consav import quadrature

import time

def check_participation_constraints(power_idx,power,Sw,Sm,idx_single,idx_couple,list_start_as_couple,list_remain_couple,list_trans_to_single, par):

    # step 1: check end points
    # 1a: check if all values are consistent with marriage
    if (Sw[0] > 0.0) & (Sm[-1] > 0.0): # all values are consistent with marriage
        remain(power_idx, power,idx_single,idx_couple, list_start_as_couple,list_remain_couple, par)

    # 1b: check if all values are consistent with divorce
    elif (Sw[-1] < 0.0) | (Sm[0] < 0.0): # no value is consistent with marriage
        divorce(power_idx, power,idx_single,idx_couple, list_start_as_couple,list_trans_to_single, par)

    # 1c: check if husband is always happy, wife has indifference point
    elif (Sw[0]<0.0) & (Sw[-1] > 0.0) & (Sm[-1] > 0.0): # husband is always happy, wife has indifference point
        # find wife's indifference point
        Low_w = find_left_point(Sw, par) + 1 # lowest point where she has positive surplus
        power_at_zero_w = linear_interp_1d._interp_1d(Sw,par.grid_power,0.0,Low_w-1) # interpolated power at indifference point

        # update case 1c
        for iP in range(par.num_power):
            idx = idx_couple(iP)

            if iP < Low_w: # update to wife's indifference point
                power_idx[idx] = Low_w
                power[idx] = power_at_zero_w
            else:
                power_idx[idx] = iP
                power[idx] = par.grid_power[iP]

            for i,key in enumerate(list_start_as_couple):
                if iP == 0:
                    list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i],power_at_zero_w,Low_w-1)
                elif iP < Low_w: # update to wife's indifference point
                    list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple(0)]    
                else: # no bargaining
                    list_start_as_couple[i][idx] = list_remain_couple[i][iP]                


    # 1d: check if wife is always happy, husband has indifference point
    elif (Sm[0]>0.0) & (Sm[-1] < 0.0) & (Sw[0] > 0.0): # wife is always happy, husband has indifference point
        # find husband's indifference point
        Low_m = find_left_point(Sm, par) # lowest point where he has positive surplus
        power_at_zero_m = linear_interp_1d._interp_1d(Sm,par.grid_power,0.0,Low_m) # interpolated power at indifference point

        # update case 1d
        for iP in range(par.num_power):
            idx = idx_couple(iP)

            if iP > Low_m:
                power_idx[idx] = Low_m
                power[idx] = power_at_zero_m
            else:
                power_idx[idx] = iP
                power[idx] = par.grid_power[iP]
            
            for i,key in enumerate(list_start_as_couple):
                if iP <= Low_m: # no bargaining
                    list_start_as_couple[i][idx] = list_remain_couple[i][iP] # no bargaining
                elif iP == Low_m+1: # update to husbands indifference point
                    list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i],power_at_zero_m,Low_m)
                else:
                    list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple(Low_m+1)]


    # 1e: Both have indifference points
    else:
        # find indifference points
        Low_w = find_left_point(Sw, par) + 1 # lowest point where she has positive surplus
        power_at_zero_w = linear_interp_1d._interp_1d(Sw,par.grid_power,0.0,Low_w-1) # interpolated power at indifference point

        Low_m = find_left_point(Sm, par) # lowest point where he has positive surplus
        power_at_zero_m = linear_interp_1d._interp_1d(Sm,par.grid_power,0.0,Low_m) # interpolated power at indifference point

        # step 3: update bargaining power and consumption allocations
        if power_at_zero_w<=power_at_zero_m: # no divorce

            # loop through range of bargaining power
            for iP in range(par.num_power): 
                idx = idx_couple(iP)

                # update bargaining power
                if iP<Low_w: # update to woman's indifference point
                    power_idx[idx] = Low_w
                    power[idx] = power_at_zero_w
                elif iP>Low_m: # update to man's indifference point
                    power_idx[idx] = Low_m
                    power[idx] = power_at_zero_m
                else: # bargaining power constant between Low_w and Low_m
                    power_idx[idx] = iP
                    power[idx] = par.grid_power[iP]

                # update consumption allocations
                for i,key in enumerate(list_start_as_couple):
                    if (iP==0): # update to woman's indifference point
                        list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i],power_at_zero_w,Low_w-1)
                    elif iP < Low_w: # bargaining power constant until Low_w
                        list_start_as_couple[i][idx]=list_start_as_couple[i][idx_couple(0)]
                    elif (iP >= Low_w) & (iP <= Low_m): # No change between Low_w and Low_m   
                        list_start_as_couple[i][idx] = list_remain_couple[i][iP]
                    elif iP == Low_m+1: # update to man's indifference point
                        list_start_as_couple[i][idx] = linear_interp_1d._interp_1d(par.grid_power,list_remain_couple[i],power_at_zero_m,Low_m)
                    else:   # bargaining power constant after Low_m
                        list_start_as_couple[i][idx] = list_start_as_couple[i][idx_couple(Low_m+1)]; # re-use that the interpolated values are identical

        else: # divorce
            divorce(power_idx, power,idx_single,idx_couple, list_start_as_couple,list_trans_to_single, par)




def divorce(power_idx, power,idx_single,idx_couple, list_start_as_couple,list_trans_to_single, par):
    for iP in range(par.num_power):
        # overwrite output for couple
        idx = idx_couple(iP)
        for i,key in enumerate(list_start_as_couple):
            list_start_as_couple[i][idx] = list_trans_to_single[i][idx_single]

        power_idx[idx] = -1
        power[idx] = -1.0


def remain(power_idx, power,idx_single,idx_couple, list_start_as_couple,list_remain_couple, par):
    for iP in range(par.num_power):
        # overwrite output for couple
        idx = idx_couple(iP)
        for i,key in enumerate(list_start_as_couple):
            list_start_as_couple[i][idx] = list_remain_couple[i][iP]

        power_idx[idx] = iP
        power[idx] = par.grid_power[iP]




def find_indifference_point(S,par,gender):

    # check that there is an indifference point
    if S[0]*S[-1] < 0:
        for iP in range(par.num_power-1):
            if S[iP]*S[iP+1] < 0:
                denom = (par.grid_power[iP+1] - par.grid_power[iP])
                ratio = (S[iP+1] - S[iP])/denom
                power_at_zero = par.grid_power[iP] - S[iP]/ratio  # interpolated power at indifference point
                return iP, power_at_zero # return point below indifference point and power at indifference point
    else:
        if gender == man:
            return par.num_power-1, par.grid_power[-1] # check
        if gender == woman:
            return 0, par.grid_power[0] # check
        

def find_left_point(S, par):
    # flip if descending
    if S[0] > S[-1]:
        S = -S
    return linear_interp.binary_search(0,len(S), S, 0.0)

