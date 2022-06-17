import numpy as np
import numba as nb
import scipy.optimize as optimize

from EconModel import EconModelClass
from consav.grids import nonlinspace
from consav import linear_interp, linear_interp_1d
from consav import quadrature

# set gender indication as globals
woman = 1
man = 2

class HouseholdModelClass(EconModelClass):
    
    def settings(self):
        """ fundamental settings """

        # a. namespaces
        self.namespaces = []
        
        # b. other attributes
        self.other_attrs = []
        
        # c. savefolder
        self.savefolder = 'saved'
        
        # d. cpp
        self.cpp_filename = 'cppfuncs/solve.cpp'
        self.cpp_options = {'compiler':'vs'}
        
    def setup(self):
        par = self.par
        
        par.R = 1.03
        par.beta = 1.0/par.R # Discount factor
        
        par.div_A_share = 0.5 # divorce share of wealth to wife
        
        # income
        par.inc_w = 1.0
        par.inc_m = 1.0

        # Utility: gender-specific parameters
        par.rho_w = 2.0        # CRRA
        par.rho_m = 2.0        # CRRA
        
        par.alpha1_w = 1.0
        par.alpha1_m = 1.0
        
        par.alpha2_w = 1.0
        par.alpha2_m = 1.0
        
        par.phi_w = 0.2
        par.phi_m = 0.2
        
        # state variables
        par.T = 10
        
        # wealth
        par.num_A = 50
        par.max_A = 10.0
        
        # bargaining power
        par.num_power = 21

        # love/match quality
        par.num_love = 11
        par.max_love = 2.0

        par.sigma_love = 0.1
        par.num_shock_love = 5

        # pre-computation
        par.num_Ctot = 50
        par.max_Ctot = par.max_A*5

        # simulation
        par.simT = par.T
        par.simN = 1000

        # cpp
        par.do_cpp = False
        par.threads = 8
        
    def allocate(self):
        par = self.par
        sol = self.sol
        sim = self.sim

        # setup grids
        self.setup_grids()
        
        # singles
        shape_single = (par.T,par.num_A)
        sol.Vw_single = np.nan + np.ones(shape_single)
        sol.Vm_single = np.nan + np.ones(shape_single)
        
        sol.Cw_priv_single = np.nan + np.ones(shape_single)
        sol.Cm_priv_single = np.nan + np.ones(shape_single)
        sol.Cw_pub_single = np.nan + np.ones(shape_single)
        sol.Cm_pub_single = np.nan + np.ones(shape_single)

        # couples
        shape_couple = (par.T,par.num_power,par.num_love,par.num_A)
        sol.Vw_couple = np.nan + np.ones(shape_couple)
        sol.Vm_couple = np.nan + np.ones(shape_couple)
        
        sol.Cw_priv_couple = np.nan + np.ones(shape_couple)
        sol.Cm_priv_couple = np.nan + np.ones(shape_couple)
        sol.C_pub_couple = np.nan + np.ones(shape_couple)

        sol.power_idx = np.zeros(shape_couple,dtype=np.int_)

        # pre-compute optimal consumption allocation
        shape_pre = (par.num_power,par.num_Ctot)
        sol.pre_Ctot_Cw_priv = np.nan + np.ones(shape_pre)
        sol.pre_Ctot_Cm_priv = np.nan + np.ones(shape_pre)
        sol.pre_Ctot_C_pub = np.nan + np.ones(shape_pre)
        
        # simulation
        shape_sim = (par.simN,par.simT)
        sim.Cw_priv = np.nan + np.ones(shape_sim)
        sim.Cm_priv = np.nan + np.ones(shape_sim)
        sim.Cw_pub = np.nan + np.ones(shape_sim)
        sim.Cm_pub = np.nan + np.ones(shape_sim)
        
        sim.A = np.nan + np.ones(shape_sim)
        sim.Aw = np.nan + np.ones(shape_sim)
        sim.Am = np.nan + np.ones(shape_sim)
        sim.couple = np.nan + np.ones(shape_sim)
        sim.power_idx = np.ones(shape_sim,dtype=np.int_)
        sim.power = np.nan + np.ones(shape_sim)
        sim.love = np.nan + np.ones(shape_sim)

        # shocks
        np.random.seed(9210)
        sim.draw_love = np.random.normal(size=shape_sim)

        # initial distribution
        sim.init_A = np.zeros(par.simN)
        sim.init_Aw = np.zeros(par.simN)
        sim.init_Am = np.zeros(par.simN)
        sim.init_couple = np.ones(par.simN,dtype=np.bool)
        sim.init_power_idx = par.num_power//2 * np.ones(par.simN,dtype=np.int_)
        sim.init_love = np.zeros(par.simN)
        
    def setup_grids(self):
        par = self.par
        
        # wealth. Single grids are such to avoid interpolation
        par.grid_A = nonlinspace(1.0e-6,par.max_A,par.num_A,1.1)

        par.grid_A_single = np.ones((2,par.num_A))
        par.grid_A_single[0] = par.div_A_share * par.grid_A
        par.grid_A_single[1] = (1.0 - par.div_A_share) * par.grid_A

        # power. non-linear grid with more mass in both tails.
        odd_num = np.mod(par.num_power,2)
        first_part = nonlinspace(0.0,0.5,(par.num_power+odd_num)//2,1.3)
        last_part = np.flip(1.0 - nonlinspace(0.0,0.5,(par.num_power-odd_num)//2 + 1,1.3))[1:]
        par.grid_power = np.append(first_part,last_part)

        # love grid and shock
        if par.num_love>1:
            par.grid_love = np.linspace(-par.max_love,par.max_love,par.num_love)
        else:
            par.grid_love = np.array([0.0])

        if par.sigma_love<=1.0e-6:
            par.num_shock_love = 1
            par.grid_shock_love,par.grid_weight_love = np.array([0.0]),np.array([1.0])
        else:
            par.grid_shock_love,par.grid_weight_love = quadrature.normal_gauss_hermite(par.sigma_love,par.num_shock_love)

        # pre-computation
        par.grid_Ctot = nonlinspace(1.0e-6,par.max_Ctot,par.num_Ctot,1.0)

    def solve(self):
        sol = self.sol
        par = self.par 

        # setup grids
        self.setup_grids()
        
        if par.do_cpp:
            self.cpp.solve(sol,par)

        else:
            
            # precompute the optimal intra-temporal consumption allocation for couples given total consumpotion
            for iP,power in enumerate(par.grid_power):
                for i,C_tot in enumerate(par.grid_Ctot):
                    sol.pre_Ctot_Cw_priv[iP,i], sol.pre_Ctot_Cm_priv[iP,i], sol.pre_Ctot_C_pub[iP,i] = solve_intraperiod_couple(C_tot,power,par)

            # loop backwards
            for t in reversed(range(par.T)):
                self.solve_single(t)
                self.solve_couple(t)

    def simulate(self):
        sol = self.sol
        sim = self.sim
        par = self.par

        for i in range(par.simN):
            for t in range(par.simT):

                # state variables
                if t==0:
                    A_lag = sim.init_A[i]
                    Aw_lag = sim.init_Aw[i]
                    Am_lag = sim.init_Am[i]
                    couple_lag = sim.init_couple[i]
                    power_idx_lag = sim.init_power_idx[i]
                    love_lag = sim.init_love[i]

                else:
                    A_lag = sim.A[i,t-1]
                    Aw_lag = sim.Aw[i,t-1]
                    Am_lag = sim.Am[i,t-1]
                    couple_lag = sim.couple[i,t-1]
                    power_idx_lag = sim.power_idx[i,t-1]
                    love_lag = sim.love[i,t-1]
                
                power_lag = par.grid_power[power_idx_lag]

                # first check if they want to remain together and what the bargaining power will be if they do.
                if couple_lag:                   
                    # love-shock
                    sim.love[i,t] = love = love_lag + par.sigma_love*sim.draw_love[i,t]

                    # use the power index:
                    power_idx = np.round(linear_interp.interp_2d(par.grid_love,par.grid_A,np.float_(sol.power_idx[t,power_idx_lag]),love,A_lag))

                    if power_idx < 0.0: # divorce is coded as -1
                        sim.couple[i,t] = False

                    else:
                        sim.couple[i,t] = True

                else: # remain single

                    sim.couple[i,t] = False


                # update behavior
                if sim.couple[i,t]:
                    
                    # optimal consumption allocation if couple
                    sol_C_tot = sol.Cw_priv_couple[t,power_idx_lag] + sol.Cm_priv_couple[t,power_idx_lag] + sol.C_pub_couple[t,power_idx_lag]
                    C_tot = linear_interp.interp_2d(par.grid_love,par.grid_A,sol_C_tot,love,A_lag)

                    sim.Cw_priv[i,t], sim.Cm_priv[i,t], C_pub = intraperiod_allocation(C_tot,power_idx_lag,sol,par)
                    sim.Cw_pub[i,t] = C_pub
                    sim.Cm_pub[i,t] = C_pub

                    # update end-of-period states
                    M_resources = par.R*A_lag + par.inc_w + par.inc_m
                    sim.A[i,t] = M_resources - sim.Cw_priv[i,t] - sim.Cm_priv[i,t] - sim.Cw_pub[i,t]

                    # in case of divorce
                    sim.Aw[i,t] = par.div_A_share * sim.A[i,t]
                    sim.Am[i,t] = (1.0-par.div_A_share) * sim.A[i,t]

                    sim.power_idx[i,t] = power_idx
                    sim.power[i,t] = par.grid_power[sim.power_idx[i,t]]

                else: # single

                    # optimal consumption allocations
                    sol_Cw_tot = sol.Cw_priv_single[t] + sol.Cw_pub_single[t]
                    sol_Cm_tot = sol.Cm_priv_single[t] + sol.Cm_pub_single[t]

                    Cw_tot = linear_interp.interp_1d(par.grid_A_single[woman-1],sol_Cw_tot,Aw_lag)
                    Cm_tot = linear_interp.interp_1d(par.grid_A_single[man-1],sol_Cm_tot,Am_lag)

                    sim.Cw_priv[i,t] = cons_priv_single(Cw_tot,woman,par)
                    sim.Cw_pub[i,t] = Cw_tot - sim.Cw_priv[i,t]
                    
                    sim.Cm_priv[i,t] = cons_priv_single(Cm_tot,man,par)
                    sim.Cm_pub[i,t] = Cm_tot - sim.Cm_priv[i,t]

                    # update end-of-period states
                    Mw = par.R*Aw_lag + par.inc_w
                    Mm = par.R*Am_lag + par.inc_m
                    sim.Aw[i,t] = Mw - sim.Cw_priv[i,t] - sim.Cw_pub[i,t]
                    sim.Am[i,t] = Mm - sim.Cm_priv[i,t] - sim.Cm_pub[i,t]

                    sim.power_idx[i,t] = -1
                    sim.power[i,t] = np.nan

                    sim.love[i,t] = np.nan
                    sim.A[i,t] = np.nan

               
    def solve_single(self,t):
        par = self.par
        sol = self.sol
        
        # terminal period
        if t == (par.T-1):
            for iA in range(par.num_A):
                idx = (t,iA)

                Aw = par.grid_A_single[0,iA]
                Am = par.grid_A_single[1,iA]

                Cw = par.R*Aw + par.inc_w
                Cm = par.R*Am + par.inc_m
                
                sol.Cw_priv_single[idx] = cons_priv_single(Cw,woman,par)
                sol.Cw_pub_single[idx] = Cw - sol.Cw_priv_single[idx]
                sol.Vw_single[idx] = util(sol.Cw_priv_single[idx],sol.Cw_pub_single[idx],woman,par)
                
                sol.Cm_priv_single[idx] = cons_priv_single(Cm,man,par)
                sol.Cm_pub_single[idx] = Cm - sol.Cm_priv_single[idx]
                sol.Vm_single[idx] = util(sol.Cm_priv_single[idx],sol.Cm_pub_single[idx],man,par)
            
        # earlier periods
        else:
            for iA in range(par.num_A):
                idx = (t,iA)

                Aw = par.grid_A_single[0,iA]
                Am = par.grid_A_single[1,iA]
                
                # resources
                Mw = par.R*Aw + par.inc_w
                Mm = par.R*Am + par.inc_m
                
                # search over optimal total consumption, C
                obj_w = lambda C_tot: - self.value_of_choice_single(C_tot[0],Mw,woman,sol.Vw_single[t+1])
                obj_m = lambda C_tot: - self.value_of_choice_single(C_tot[0],Mm,man,sol.Vm_single[t+1])
                
                res_w = optimize.minimize(obj_w,Mw/2.0,bounds=((1.0e-8,Mw),))
                res_m = optimize.minimize(obj_m,Mm/2.0,bounds=((1.0e-8,Mm),))
                
                # store results
                Cw = res_w.x
                sol.Cw_priv_single[idx] = cons_priv_single(Cw,woman,par)
                sol.Cw_pub_single[idx] = Cw - sol.Cw_priv_single[idx]
                
                Cm = res_m.x
                sol.Cm_priv_single[idx] = cons_priv_single(Cm,man,par)
                sol.Cm_pub_single[idx] = Cm - sol.Cm_priv_single[idx]
                
                sol.Vw_single[idx] = -res_w.fun
                sol.Vm_single[idx] = -res_m.fun

    def solve_couple(self,t):
        par = self.par
        sol = self.sol

        tmp_Vw,tmp_Vm,tmp_Cw_priv,tmp_Cm_priv,tmp_C_pub = np.ones(par.num_power),np.ones(par.num_power),np.ones(par.num_power),np.ones(par.num_power),np.ones(par.num_power)
        
        Vw_next = None
        Vm_next = None
        for iL,love in enumerate(par.grid_love):
            for iA,A in enumerate(par.grid_A):
                M_resources = par.R*A + par.inc_w + par.inc_m
                
                starting_val = None
                for iP,power in enumerate(par.grid_power):
                    # continuation values
                    if t<(par.T-1):
                        Vw_next = self.sol.Vw_couple[t+1,iP]
                        Vm_next = self.sol.Vm_couple[t+1,iP]

                    # starting values
                    if iP>0:
                        C_tot_last = tmp_Cw_priv[iP-1] + tmp_Cm_priv[iP-1] + tmp_C_pub[iP-1]
                        starting_val = np.array([C_tot_last])
                    
                    # solve unconstrained problem
                    tmp_Cw_priv[iP], tmp_Cm_priv[iP], tmp_C_pub[iP], tmp_Vw[iP], tmp_Vm[iP] = self.solve_uncon_couple(t,M_resources,iL,iP,power,Vw_next,Vm_next,starting_val=starting_val)

                # check the participation constraints
                idx_single = (t,iA)
                idx_couple = lambda iP: (t,iP,iL,iA)

                list_couple = (sol.Vw_couple,sol.Vm_couple,sol.Cw_priv_couple,sol.Cm_priv_couple,sol.C_pub_couple)
                list_raw = (tmp_Vw,tmp_Vm,tmp_Cw_priv,tmp_Cm_priv,tmp_C_pub)
                list_single = (sol.Vw_single,sol.Vm_single,sol.Cw_priv_single,sol.Cm_priv_single,sol.Cw_pub_single) # last input here not important in case of divorce
                
                Sw = tmp_Vw - sol.Vw_single[idx_single] 
                Sm = tmp_Vm - sol.Vm_single[idx_single] 
                
                check_participation_constraints(sol.power_idx,Sw,Sm,idx_single,idx_couple,list_couple,list_raw,list_single, par)

    def solve_uncon_couple(self,t,M_resources,iL,iP,power,Vw_next,Vm_next,starting_val = None):
        par = self.par

        if t==(par.T-1): # Terminal period
            C_tot = M_resources

        else:
            # objective function
            obj = lambda x: - self.value_of_choice_couple(x[0],t,M_resources,iL,iP,power,Vw_next,Vm_next)[0]

            # bound: credit constraint, no borrowing
            bounds = optimize.Bounds(1.0e-6, M_resources - 1.0e-6, keep_feasible=True)

            # optimize
            x0 = np.array([M_resources * 0.8]) if starting_val is None else starting_val
            res = optimize.minimize(obj,x0,bounds=bounds)
            C_tot = res.x[0]

        # implied consumption allocation (re-calculation)
        _, Cw_priv, Cm_priv, C_pub, Vw,Vm = self.value_of_choice_couple(C_tot,t,M_resources,iL,iP,power,Vw_next,Vm_next)

        # return objects
        return Cw_priv, Cm_priv, C_pub, Vw, Vm

    def value_of_choice_couple(self,C_tot,t,M_resources,iL,iP,power,Vw_next,Vm_next):
        sol = self.sol
        par = self.par

        love = par.grid_love[iL]
        
        # current utility from consumption allocation
        Cw_priv, Cm_priv, C_pub = intraperiod_allocation(C_tot,iP,sol,par)
        Vw = util(Cw_priv,C_pub,woman,par,love)
        Vm = util(Cm_priv,C_pub,man,par,love)

        # add continuation value [TODO: re-use index would speed this up since only output different!]
        if t < (par.T-1):
            savings = M_resources - C_tot 
            EVw_plus = 0.0
            EVm_plus = 0.0
            for iL_next in range(par.num_shock_love):
                love_next = love + par.grid_shock_love[iL_next]

                EVw_plus += par.grid_weight_love[iL_next] * linear_interp.interp_2d(par.grid_love,par.grid_A , Vw_next, love_next,savings)
                EVm_plus += par.grid_weight_love[iL_next] * linear_interp.interp_2d(par.grid_love,par.grid_A , Vm_next, love_next,savings)

            Vw += par.beta*EVw_plus
            Vm += par.beta*EVm_plus

        # return
        Val = power*Vw + (1.0-power)*Vm
        return Val , Cw_priv, Cm_priv, C_pub, Vw,Vm


    def value_of_choice_couple_old(self,C_tot,t,M_resources,iL,iP,power,Vw_next,Vm_next):
        sol = self.sol
        par = self.par

        love = par.grid_love[iL]
        
        # current utility from consumption allocation (intra-period)
        Cw_priv, Cm_priv, C_pub = intraperiod_allocation(C_tot,iP,sol,par)

        if t==(par.T-1): # terminal period
            Vw = util(Cw_priv,C_pub,woman,par,love)
            Vm = util(Cm_priv,C_pub,man,par,love)

        else:
            # value of choice of household
            savings = M_resources - C_tot 
            Vw = indiv_value_of_choice_couple_old(Cw_priv,C_pub, woman,savings,love, Vw_next, par)
            Vm = indiv_value_of_choice_couple_old(Cm_priv,C_pub, man,savings,love, Vm_next, par)

        # return
        Val = power*Vw + (1.0-power)*Vm
        return Val , Cw_priv, Cm_priv, C_pub, Vw,Vm
        
    def value_of_choice_single(self,C_tot,M,gender,V_next):
        par = self.par

        # flow-utility
        C_priv = cons_priv_single(C_tot,gender,par)
        C_pub = C_tot - C_priv
        
        Util = util(C_priv,C_pub,gender,par)
        
        # continuation value
        grid_A = par.grid_A_single[gender-1]
        A = M - C_tot

        Vnext = linear_interp.interp_1d(grid_A,V_next,A)
        
        # return discounted sum
        return Util + par.beta*Vnext
   

def util(c_priv,c_pub,gender,par,love=0.0):
    if gender == 1:
        rho = par.rho_w
        phi = par.phi_w
        alpha1 = par.alpha1_w
        alpha2 = par.alpha2_w
    else:
        rho = par.rho_m
        phi = par.phi_m
        alpha1 = par.alpha1_m
        alpha2 = par.alpha2_m
    
    return ((alpha1*c_priv**phi + alpha2*c_pub**phi)**(1.0-rho))/(1.0-rho) + love

##########
# Single #
##########
def cons_priv_single(C_tot,gender,par):
    # closed form solution for intra-period problem of single.
    if gender == 1:
        rho = par.rho_w
        phi = par.phi_w
        alpha1 = par.alpha1_w
        alpha2 = par.alpha2_w
    else:
        rho = par.rho_m
        phi = par.phi_m
        alpha1 = par.alpha1_m
        alpha2 = par.alpha2_m   
    
    return C_tot/(1.0 + (alpha2/alpha1)**(1.0/(1.0-phi)) )

##########
# Couple #
##########   
def intraperiod_allocation(C_tot,iP,sol,par):

    # interpolate pre-computed solution (do smarter if working)
    j1 = linear_interp.binary_search(0,par.num_Ctot,par.grid_Ctot,C_tot)
    Cw_priv = linear_interp_1d._interp_1d(par.grid_Ctot,sol.pre_Ctot_Cw_priv[iP],C_tot,j1)
    Cm_priv = linear_interp_1d._interp_1d(par.grid_Ctot,sol.pre_Ctot_Cm_priv[iP],C_tot,j1)
    C_pub = linear_interp_1d._interp_1d(par.grid_Ctot,sol.pre_Ctot_C_pub[iP],C_tot,j1) 

    return Cw_priv, Cm_priv, C_pub

def solve_intraperiod_couple(C_tot,power,par,starting_val=None):
    
    # setup estimation. Impose constraint that C_tot = Cw+Cm+C
    bounds = optimize.Bounds(0.0, C_tot, keep_feasible=True)
    obj = lambda x: - (power*util(x[0],C_tot-np.sum(x),1,par) + (1.0-power)*util(x[1],C_tot-np.sum(x),2,par))
    
    # estimate
    x0 = np.array([C_tot/3,C_tot/3]) if starting_val is None else starting_val
    res = optimize.minimize(obj,x0,bounds=bounds)
    # constr = {'type':'ineq','fun':lambda x: C_tot - np.sum(x)}
    # res = optimize.minimize(obj,x0,bounds=bounds,constraints=constr)

    # unpack
    Cw_priv = res.x[0]
    Cm_priv = res.x[1]
    C_pub = C_tot - Cw_priv - Cm_priv

    return Cw_priv,Cm_priv,C_pub

def indiv_value_of_choice_couple_old(cons_priv,cons_pub, gender,savings,love, V_next, par):
    # value of choice for an individual in a couple
    # savings: household level savings
    # V_next: value of starting in a couple in the beginning of next period

    # current utility
    Util = util(cons_priv,cons_pub,gender,par,love)

    # continuation value
    EV_plus = 0.0
    for iL_next in range(par.num_shock_love):
        love_next = love + par.grid_shock_love[iL_next]
        EV_plus += par.grid_weight_love[iL_next] * linear_interp.interp_2d(par.grid_love,par.grid_A , V_next, love_next,savings)

    # discounted sum
    return Util + par.beta*EV_plus

def check_participation_constraints(power,Sw,Sm,idx_single,idx_couple,list_couple,list_raw,list_single, par):
    
    # check the participation constraints. Array
    min_Sw = np.min(Sw)
    min_Sm = np.min(Sm)
    max_Sw = np.max(Sw)
    max_Sm = np.max(Sm)

    if (min_Sw >= 0.0) & (min_Sm >= 0.0): # all values are consistent with marriage
        for iP in range(par.num_power):

            # overwrite output for couple
            idx = idx_couple(iP)
            for i,key in enumerate(list_couple):
                list_couple[i][idx] = list_raw[i][iP]

            power[idx] = iP

    elif (max_Sw < 0.0) | (max_Sm < 0.0): # no value is consistent with marriage
        for iP in range(par.num_power):

            # overwrite output for couple
            idx = idx_couple(iP)
            for i,key in enumerate(list_couple):
                list_couple[i][idx] = list_single[i][idx_single]

            power[idx] = -1

    else: 
    
        # find lowest (highest) value with positive surplus for women (men)
        Low_w = 0 # in case there is no crossing, this will be the correct value
        Low_m = par.num_power-1 # in case there is no crossing, this will be the correct value
        for iP in range(par.num_power-1):
            if (Sw[iP]<0) & (Sw[iP+1]>=0):
                Low_w = iP+1
                
            if (Sm[iP]>=0) & (Sm[iP+1]<0):
                Low_m = iP

        # update the outcomes
        for iP in range(par.num_power):

            # index to store solution for couple 
            idx = idx_couple(iP)
    
            # woman wants to leave
            if iP<Low_w: 
                if Sm[Low_w] > 0: # man happy to shift some bargaining power

                    for i,key in enumerate(list_couple):
                        list_couple[i][idx] = list_raw[i][Low_w]

                    power[idx] = Low_w
                    
                else: # divorce

                    for i,key in enumerate(list_couple):
                        list_couple[i][idx] = list_single[i][idx_single]

                    power[idx] = -1
                
            # man wants to leave
            elif iP>Low_m: 
                if Sw[Low_m] > 0: # woman happy to shift some bargaining power
                    
                    for i,key in enumerate(list_couple):
                        list_couple[i][idx] = list_raw[i][Low_m]

                    power[idx] = Low_m
                    
                else: # divorce

                    for i,key in enumerate(list_couple):
                        list_couple[i][idx] = list_single[i][idx_single]

                    power[idx] = -1

            else: # no-one wants to leave

                for i,key in enumerate(list_couple):
                    list_couple[i][idx] = list_raw[i][iP]

                power[idx] = iP



