import numpy as np
import numba as nb
import scipy.optimize as optimize
import scipy.stats as stats

from EconModel import EconModelClass
from consav.grids import nonlinspace
from consav import linear_interp, linear_interp_1d
from consav import quadrature

import bargaining_algorithm as ba

import time

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
        par.div_cost = 0.0
        
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
        par.max_A = 15.0
        
        # bargaining power
        par.num_power = 21

        # love/match quality
        par.num_love = 41
        par.max_love = 1.0

        par.sigma_love = 0.1
        par.num_shock_love = 5

        # re-partnering
        par.p_meet = 0.1
        par.prob_partner_A_w = np.array([[np.nan]]) # if not set here, defaults to np.eye(par.num_A) in setup_grids
        par.prob_partner_A_m = np.array([[np.nan]])

        # pre-computation
        par.interp_inverse = False # True: interpolate inverse consumption

        par.num_Ctot = 100
        par.max_Ctot = par.max_A*2
        
        par.do_egm = False
        par.analytic_marg_u_single = False
        par.analytic_inv_marg_u_single = False
        par.num_A_pd = par.num_A * 2
        par.max_A_pd = par.max_A
        par.num_marg_u = 200

        # simulation
        par.seed = 9210
        par.simT = par.T
        par.simN = 50_000
                
        # cpp
        par.do_cpp = False
        par.threads = 8

        par.interp_power = False
        par.centered_gradient = False
        
    def allocate(self):
        par = self.par
        sol = self.sol
        sim = self.sim

        # setup grids
        self.setup_grids()
        
        # a. singles
        shape_single = (par.T,par.num_A)                        # single states: T and assets

        # a.1. single to single
        sol.Vw_single_to_single = np.nan + np.ones(shape_single)
        sol.Vm_single_to_single = np.nan + np.ones(shape_single) 

        sol.Cw_priv_single_to_single = np.nan + np.ones(shape_single)     # private consumption, single
        sol.Cm_priv_single_to_single = np.nan + np.ones(shape_single)
        sol.Cw_pub_single_to_single = np.nan + np.ones(shape_single)      # public consumption, single (i guess this functions as a private good when single)
        sol.Cm_pub_single_to_single = np.nan + np.ones(shape_single)
        sol.Cw_tot_single_to_single = np.nan + np.ones(shape_single)      # total concumption, single
        sol.Cm_tot_single_to_single = np.nan + np.ones(shape_single)

        ### a.1.1. post-decision grids (EGM)
        sol.EmargUw_single_to_single_pd = np.zeros(shape_single)           # Expected marginal utility post-decision, woman single
        sol.C_totw_single_to_single_pd = np.zeros(par.num_A_pd)            # C for EGM, woman single 
        sol.Mw_single_to_single_pd = np.zeros(par.num_A_pd)                # Endogenous grid, woman single

        sol.EmargUm_single_to_single_pd = np.zeros(shape_single)          # Expected marginal utility post-decision, man single
        sol.C_totm_single_to_single_pd = np.zeros(par.num_A_pd)           # C for EGM, man single
        sol.Mm_single_to_single_pd = np.zeros(par.num_A_pd)               # Endogenous grid, man single

        ## a.2. couple to single
        sol.Vw_couple_to_single = np.nan + np.ones(shape_single)        # Value marriage -> single
        sol.Vm_couple_to_single = np.nan + np.ones(shape_single)

        sol.Cw_priv_couple_to_single = np.nan + np.ones(shape_single)   # Private consumption marriage -> single
        sol.Cm_priv_couple_to_single = np.nan + np.ones(shape_single)
        sol.Cw_pub_couple_to_single = np.nan + np.ones(shape_single)    # Public consumption marriage -> single 
        sol.Cm_pub_couple_to_single = np.nan + np.ones(shape_single)    # Not used
        sol.Cw_tot_couple_to_single = np.nan + np.ones(shape_single)
        sol.Cm_tot_couple_to_single = np.nan + np.ones(shape_single)

        ## a.3. start as single
        sol.EVw_start_as_single = np.nan + np.ones(shape_single)
        sol.EVm_start_as_single = np.nan + np.ones(shape_single)  
        sol.EmargVw_start_as_single = np.nan + np.ones(shape_single)
        sol.EmargVm_start_as_single = np.nan + np.ones(shape_single)  


        # b. couples
        shape_couple = (par.T,par.num_power,par.num_love,par.num_A)     # states when couple: T, assets, power, love

        ## b.1. couple to couple
        sol.Vw_couple_to_couple = np.nan + np.ones(shape_couple)           # value marriage -> marriage
        sol.Vm_couple_to_couple = np.nan + np.ones(shape_couple)
        sol.V_couple_to_couple = -np.inf + np.ones(shape_couple)           # value marriage -> marriage, weighted by bargaining power (for DC-EGM)
                                                                           # AMO: initialized at -inf which is useful in upper envelope
        sol.Cw_priv_couple_to_couple = np.nan + np.ones(shape_couple)      # private consumption, marriage -> marriage
        sol.Cm_priv_couple_to_couple = np.nan + np.ones(shape_couple)      
        sol.C_pub_couple_to_couple = np.nan + np.ones(shape_couple)        # public consumption, marriage -> marriage
        sol.C_tot_couple_to_couple = np.nan + np.ones(shape_couple)        # total consumption, marriage -> marriage
        
        sol.Sw = np.nan + np.ones(par.num_power)                         # Surplus of mariage
        sol.Sm = np.nan + np.ones(par.num_power)

        sol.power_idx = np.zeros(shape_couple,dtype=np.int_)            # index of bargaining weight (approx)
        sol.power = np.zeros(shape_couple)                              # bargaining weight (interpolated)

        ### b.1.1. post-decision grids (EGM)
        shape_egm = (par.T, par.num_power,par.num_love,par.num_A_pd)
        sol.EmargU_pd = np.zeros(shape_egm)                     # Expected marginal utility post-decision
        sol.C_tot_pd = np.zeros(shape_egm)                      # C for EGM
        sol.M_pd = np.zeros(shape_egm)                          # Endogenous grid
        sol.V_couple_to_couple_pd = np.zeros(shape_egm)         # Value of being couple, post-decision

        ## b.2. single to couple
        sol.Vw_single_to_couple = np.nan + np.ones(shape_couple)           # value single -> marriage
        sol.Vm_single_to_couple = np.nan + np.ones(shape_couple)
        sol.V_single_to_couple = -np.inf + np.ones(shape_couple)           
                                                                        
        sol.Cw_priv_single_to_couple = np.nan + np.ones(shape_couple)      
        sol.Cm_priv_single_to_couple = np.nan + np.ones(shape_couple)      
        sol.C_pub_single_to_couple = np.nan + np.ones(shape_couple)        
        sol.Cw_tot_single_to_couple = np.nan + np.ones(shape_couple)   
        sol.Cm_tot_single_to_couple = np.nan + np.ones(shape_couple) 
  
        shape_power =(par.T,par.num_love,par.num_A,par.num_A)          
        sol.initial_power = np.nan + np.zeros(shape_power)
        sol.initial_power_idx = np.zeros(shape_power,dtype=np.int_)

        ## b.3. start as couple
        sol.Vw_start_as_couple = np.nan + np.ones(shape_couple)
        sol.Vm_start_as_couple = np.nan + np.ones(shape_couple)
        sol.margV_start_as_couple = np.zeros(shape_couple) 

        sol.Cw_priv_start_as_couple = np.nan + np.ones(shape_couple)             # private consumption, couple
        sol.Cm_priv_start_as_couple = np.nan + np.ones(shape_couple)             
        sol.C_pub_start_as_couple = np.nan + np.ones(shape_couple)               # public consumption, couple
        sol.C_tot_start_as_couple = np.zeros(shape_couple)                       # total consumption, couple
                                                                        # AMO: initialized at zero which is useful in upper envelope
        
        # c. temporary containers
        sol.savings_vec = np.zeros(par.num_shock_love)          
        sol.Vw_plus_vec = np.zeros(par.num_shock_love)          # not sure (next period maybe)
        sol.Vm_plus_vec = np.zeros(par.num_shock_love) 


        # d. pre-compute optimal consumption allocation
        shape_pre = (par.num_power,par.num_Ctot)
        sol.pre_Ctot_Cw_priv = np.nan + np.ones(shape_pre)      # precomputed optimal allocation of consumption over grid of total C 
        sol.pre_Ctot_Cm_priv = np.nan + np.ones(shape_pre)
        sol.pre_Ctot_C_pub = np.nan + np.ones(shape_pre)


        # e. simulation
        # NB: all arrays not containing "init" or "draw" in name are wiped before each simulation
        shape_sim = (par.simN,par.simT)
        sim.Cw_priv = np.nan + np.ones(shape_sim)               
        sim.Cm_priv = np.nan + np.ones(shape_sim)
        sim.Cw_pub = np.nan + np.ones(shape_sim)
        sim.Cm_pub = np.nan + np.ones(shape_sim)
        sim.Cw_tot = np.nan + np.ones(shape_sim)
        sim.Cm_tot = np.nan + np.ones(shape_sim)
        sim.C_tot = np.nan + np.ones(shape_sim)
        
        sim.A = np.nan + np.ones(shape_sim)
        sim.Aw = np.nan + np.ones(shape_sim)
        sim.Am = np.nan + np.ones(shape_sim)
        sim.couple = np.nan + np.ones(shape_sim)
        sim.power = np.nan + np.ones(shape_sim)
        sim.love = np.nan + np.ones(shape_sim)

        # containers for verifying simulaton
        sim.A_own = np.nan + np.ones(shape_sim)
        sim.A_partner = np.nan + np.ones(shape_sim)

        ## e.1. shocks
        np.random.seed(par.seed)
        sim.draw_love = np.random.normal(size=shape_sim)
        sim.draw_meet = np.random.uniform(size=shape_sim) # for meeting a partner

        sim.draw_uniform_partner_Aw = np.random.uniform(size=shape_sim) # for inverse cdf transformation of partner wealth
        sim.draw_uniform_partner_Am = np.random.uniform(size=shape_sim) # for inverse cdf tranformation of partner wealth

        sim.draw_repartner_iL = np.random.choice(par.num_love, p=par.prob_partner_love, size=shape_sim) # Love index when repartnering

        ## e.2. initial distribution
        sim.init_A = par.grid_A[10] + np.zeros(par.simN)
        sim.init_Aw = sim.init_A * par.div_A_share
        sim.init_Am = sim.init_A * (1.0 - par.div_A_share)
        sim.init_couple = np.ones(par.simN,dtype=np.bool_)
        sim.init_power_idx = par.num_power//2 * np.ones(par.simN,dtype=np.int_)
        sim.init_love = np.zeros(par.simN)
        
        # f. timing
        sol.solution_time = np.array([0.0])
        
    def setup_grids(self):
        par = self.par
        
        # wealth. Single grids are such to avoid interpolation
        par.grid_A = nonlinspace(0.0,par.max_A,par.num_A,1.1)       # asset grid

        par.grid_Aw = par.div_A_share * par.grid_A                  # asset grid in case of divorce
        par.grid_Am = (1.0 - par.div_A_share) * par.grid_A

        # power. non-linear grid with more mass in both tails.
        odd_num = np.mod(par.num_power,2)
        first_part = nonlinspace(0.0,0.5,(par.num_power+odd_num)//2,1.3)
        last_part = np.flip(1.0 - nonlinspace(0.0,0.5,(par.num_power-odd_num)//2 + 1,1.3))[1:]
        par.grid_power = np.append(first_part,last_part)
        par.grid_power_flip = np.flip(par.grid_power) # flip for men

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
        par.grid_Ctot = nonlinspace(1.0e-6,par.max_Ctot,par.num_Ctot,1.1)   

        # EGM
        par.grid_util = np.nan + np.ones((par.num_power,par.num_marg_u))
        par.grid_marg_u = np.nan + np.ones(par.grid_util.shape)
        par.grid_marg_u_for_inv = np.nan + np.ones(par.grid_util.shape)

        par.grid_C_for_marg_u = nonlinspace(1.0e-6,par.max_Ctot,par.num_marg_u,1.1)

        par.grid_inv_marg_u = np.flip(par.grid_C_for_marg_u) # Flipped to make interpolation possible ## AMO: invert
        if par.interp_inverse:
            par.grid_inv_marg_u = 1.0/par.grid_inv_marg_u

        par.grid_marg_u_single_w = np.nan + np.ones((par.num_marg_u))
        par.grid_marg_u_single_w_for_inv = np.nan + np.ones((par.num_marg_u))

        par.grid_marg_u_single_m = np.nan + np.ones((par.num_marg_u))
        par.grid_marg_u_single_m_for_inv = np.nan + np.ones((par.num_marg_u))

        par.grid_A_pd = nonlinspace(0.0,par.max_A_pd,par.num_A_pd,1.1)
        par.grid_Aw_pd = par.div_A_share * par.grid_A_pd
        par.grid_Am_pd = (1.0 - par.div_A_share) * par.grid_A_pd

        # re-partering probabilities
        par.prob_repartner = par.p_meet*np.ones(par.T) # likelihood of meeting a partner

        if np.isnan(par.prob_partner_A_w[0,0]):
            par.prob_partner_A_w = np.eye(par.num_A) #np.ones((par.num_A,par.num_A))/par.num_A # likelihood of meeting a partner with a particular level of wealth, conditional on own
    
        if np.isnan(par.prob_partner_A_m[0,0]):
            par.prob_partner_A_m = np.eye(par.num_A) #np.ones((par.num_A,par.num_A))/par.num_A # likelihood of meeting a partner with a particular level of wealth, conditional on own

        # Norm distributed initial love - note: Probability mass between points (approximation of continuous distribution)
        if par.sigma_love<=1.0e-6:
            love_cdf = np.where(par.grid_love>=0.0,1.0,0.0)
        else:
            love_cdf = stats.norm.cdf(par.grid_love,0.0,par.sigma_love)
        par.prob_partner_love = np.diff(love_cdf,1)
        par.prob_partner_love = np.append(par.prob_partner_love,0.0) # lost last point in diff
        # par.prob_partner_love = np.ones(par.num_love)/par.num_love # uniform

        par.cdf_partner_Aw = np.cumsum(par.prob_partner_A_w,axis=1) # cumulative distribution to be used in simulation
        par.cdf_partner_Am = np.cumsum(par.prob_partner_A_m,axis=1)


    def solve(self):
        
        #Begin timing
        start_time = time.time()
        
        sol = self.sol
        par = self.par 

        # setup grids
        self.setup_grids()

        #######################
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

        
        # End timing
        sol.solution_time[0] = time.time() - start_time

    def simulate(self):
        sol = self.sol
        sim = self.sim
        par = self.par

        # clear simulation
        for key, val in sim.__dict__.items():
            if 'init' in key or 'draw' in key: continue
            setattr(sim, key, np.zeros(val.shape)+np.nan)

        if par.do_cpp:
            self.cpp.simulate(sim,sol,par)

        else:
            simulate(sim,sol,par)

        # total consumption
        sim.Cw_tot = sim.Cw_priv + sim.Cw_pub
        sim.Cm_tot = sim.Cm_priv + sim.Cm_pub
        sim.C_tot = sim.Cw_priv + sim.Cm_priv + sim.Cw_pub

               
    def solve_single(self,t): # vfi solution conditional on being single
        par = self.par
        sol = self.sol
        
        # loop through state variable: wealth
        for iA in range(par.num_A):
            idx = (t,iA)

            # resources
            Aw = par.grid_Aw[iA]
            Am = par.grid_Am[iA]

            Mw = resources_single(Aw,woman,par) 
            Mm = resources_single(Am,man,par) 

            if t == (par.T-1): # terminal period
                
                # intra-period allocation: consume all resources
                sol.Cw_priv_single[idx],sol.Cw_pub_single[idx] = intraperiod_allocation_single(Mw,woman,par)
                sol.Vw_single[idx] = util(sol.Cw_priv_single[idx],sol.Cw_pub_single[idx],woman,par)
                
                sol.Cm_priv_single[idx],sol.Cm_pub_single[idx] = intraperiod_allocation_single(Mm,man,par)
                sol.Vm_single[idx] = util(sol.Cm_priv_single[idx],sol.Cm_pub_single[idx],man,par)
            
            else: # earlier periods

                # search over optimal total consumption, C
                obj_w = lambda C_tot: - self.value_of_choice_single(C_tot[0],Mw,woman,sol.Vw_single[t+1])
                obj_m = lambda C_tot: - self.value_of_choice_single(C_tot[0],Mm,man,sol.Vm_single[t+1])
                
                res_w = optimize.minimize(obj_w,Mw/2.0,bounds=((1.0e-8,Mw),))
                res_m = optimize.minimize(obj_m,Mm/2.0,bounds=((1.0e-8,Mm),))
                
                # store results
                Cw = res_w.x
                sol.Cw_priv_single[idx],sol.Cw_pub_single[idx] = intraperiod_allocation_single(Cw,woman,par)
                sol.Vw_single[idx] = -res_w.fun
                
                Cm = res_m.x
                sol.Cm_priv_single[idx],sol.Cm_pub_single[idx] = intraperiod_allocation_single(Cm,man,par)
                sol.Vm_single[idx] = -res_m.fun                
                
    def solve_couple(self,t): # vfi solution conditional on being married
        

        # a. unpack and allocate temporary memory
        par = self.par
        sol = self.sol

        # b. loop through state variables
        for iL,love in enumerate(par.grid_love):
            for iA,A in enumerate(par.grid_A):
                # Construct index function and lists of arrays to be used to detrmine outcomes
                idx_single = (t,iA)
                idx_couple = lambda iP: (t,iP,iL,iA)
                
                # Calculate resources
                M_resources = resources_couple(A,par) 
                
                starting_val = None # trick :)
                for iP,power in enumerate(par.grid_power):
                    
                    # i. continuation values
                    if t<(par.T-1):
                        Vw_next = self.sol.Vw_couple[t+1,iP]
                        Vm_next = self.sol.Vm_couple[t+1,iP]
                    else:
                        Vw_next = None
                        Vm_next = None

                    # ii. starting value for total consumption. The intra-temporal allocation problem solved conditional on this
                    if iP>0:
                        C_tot_last = sol.Cw_priv_remain_couple[idx_couple(iP-1)] + sol.Cm_priv_remain_couple[idx_couple(iP-1)] + sol.C_pub_remain_couple[idx_couple(iP-1)]
                        starting_val = np.array([C_tot_last])
                    
                    # iii. solve problem if remaining married, m->m
                    sol.Cw_priv_remain_couple[idx_couple(iP)], sol.Cm_priv_remain_couple[idx_couple(iP)], sol.C_pub_remain_couple[idx_couple(iP)], sol.Vw_remain_couple[idx_couple(iP)], sol.Vm_remain_couple[idx_couple(iP)] = self.solve_remain_couple(t,M_resources,iL,iP,power,Vw_next,Vm_next,starting_val=starting_val)

                    # iv. Calculate marital surplus
                    sol.Sw[iP] = sol.Vw_remain_couple[idx_couple(iP)] - sol.Vw_single[idx_single]
                    sol.Sm[iP] = sol.Vm_remain_couple[idx_couple(iP)] - sol.Vm_single[idx_single]
                    
                # c. check the participation constraints. 
                list_start_as_couple = (sol.Vw_couple, sol.Vm_couple, sol.Cw_priv_couple, sol.Cm_priv_couple, sol.C_pub_couple)
                list_remain_couple =   (sol.Vw_remain_couple, sol.Vm_remain_couple, sol.Cw_priv_remain_couple, sol.Cm_priv_remain_couple, sol.C_pub_remain_couple)
                list_trans_to_single = (sol.Vw_single[idx_single], sol.Vm_single[idx_single], sol.Cw_priv_single[idx_single], sol.Cm_priv_single[idx_single], sol.Cw_pub_single[idx_single]) # last input here not important in case of divorce
                
                # d. check participation constraints and update sol values
                ba.check_participation_constraints(sol.power_idx,sol.power,sol.Sw,sol.Sm,idx_couple,list_start_as_couple,list_remain_couple,list_trans_to_single, par)

    def solve_remain_couple(self,t,M_resources,iL,iP,power,Vw_next,Vm_next,starting_val = None):
        par = self.par

        # a. determine optimal total consumption
        if t==(par.T-1): # Terminal period: consume all resources
            C_tot = M_resources

        else:
            # i. objective function
            obj = lambda x: - self.value_of_choice_couple(x[0],t,M_resources,iL,iP,power,Vw_next,Vm_next)[0]
            x0 = np.array([M_resources * 0.8]) if starting_val is None else starting_val

            # ii. optimize wrt. total consumption
            res = optimize.minimize(obj,x0,bounds=((1.0e-6, M_resources - 1.0e-6),) ,method='SLSQP') 
            C_tot = res.x[0]

        # b. implied intra-temporal consumption allocation (re-calculation)
        _, Cw_priv, Cm_priv, C_pub, Vw,Vm = self.value_of_choice_couple(C_tot,t,M_resources,iL,iP,power,Vw_next,Vm_next)

        # c. return objects
        return Cw_priv, Cm_priv, C_pub, Vw, Vm

    def value_of_choice_couple(self,C_tot,t,M_resources,iL,iP,power,Vw_next,Vm_next):
        # This is the value of a given total consumption choice, conditional on remaining a couple.
        sol = self.sol
        par = self.par

        love = par.grid_love[iL]
        
        # a. current utility from consumption allocation
        Cw_priv, Cm_priv, C_pub = intraperiod_allocation(C_tot,iP,sol,par)
        Vw = util(Cw_priv,C_pub,woman,par,love)
        Vm = util(Cm_priv,C_pub,man,par,love)

        # b. add continuation value if not last period
        if t < (par.T-1):
            
            # i. interpolate for vector of love shocks. Savings constant across these shocks
            sol.savings_vec[:] = M_resources - C_tot 
            love_next_vec = love + par.grid_shock_love

            linear_interp.interp_2d_vec(par.grid_love,par.grid_A , Vw_next, love_next_vec,sol.savings_vec,sol.Vw_plus_vec)
            linear_interp.interp_2d_vec(par.grid_love,par.grid_A , Vm_next, love_next_vec,sol.savings_vec,sol.Vm_plus_vec)

            # ii. expected continuation value
            EVw_plus = sol.Vw_plus_vec @ par.grid_weight_love
            EVm_plus = sol.Vm_plus_vec @ par.grid_weight_love

            # iii. add discounted expected continuation value to flow utility
            Vw += par.beta*EVw_plus
            Vm += par.beta*EVm_plus

        # c. return weighted household value and other items
        Val = power*Vw + (1.0-power)*Vm
        return Val , Cw_priv, Cm_priv, C_pub, Vw,Vm
        
    def value_of_choice_single(self,C_tot,M,gender,V_next):
        par = self.par

        # flow-utility
        C_priv = cons_priv_single(C_tot,gender,par)
        C_pub = C_tot - C_priv
        
        Util = util(C_priv,C_pub,gender,par)
        
        # continuation value
        grid_A = par.grid_Aw if gender==woman else par.grid_Am
        A = M - C_tot

        Vnext = linear_interp.interp_1d(grid_A,V_next,A)
        
        # return discounted sum
        return Util + par.beta*Vnext
   

def simulate(sim,sol,par):
    for i in range(par.simN):
        for t in range(par.simT):

            # state variables
            if t==0:
                A_lag = sim.init_A[i]
                Aw_lag = sim.init_Aw[i]
                Am_lag = sim.init_Am[i]
                couple_lag = sim.init_couple[i]
                power_idx_lag = sim.init_power_idx[i]
                love = sim.love[i,t] = sim.init_love[i]

            else:
                A_lag = sim.A[i,t-1]
                Aw_lag = sim.Aw[i,t-1]
                Am_lag = sim.Am[i,t-1]
                couple_lag = sim.couple[i,t-1]
                power_idx_lag = sim.power_idx[i,t-1]
                love = sim.love[i,t]
            
            power_lag = par.grid_power[power_idx_lag] 

            # first check if they want to remain together and what the bargaining power will be if they do.
            if couple_lag:                   

                # value of transitioning into singlehood
                Vw_single = linear_interp.interp_1d(par.grid_Aw,sol.Vw_trans_single[t],Aw_lag)
                Vm_single = linear_interp.interp_1d(par.grid_Am,sol.Vm_trans_single[t],Am_lag)

                idx = (t,power_idx_lag)
                Vw_couple_i = linear_interp.interp_2d(par.grid_love,par.grid_A,sol.Vw_remain_couple[idx],love,A_lag)
                Vm_couple_i = linear_interp.interp_2d(par.grid_love,par.grid_A,sol.Vm_remain_couple[idx],love,A_lag)

                if ((Vw_couple_i>=Vw_single) & (Vm_couple_i>=Vm_single)):
                    power_idx = power_idx_lag                            # no bargaining takes place if both have positive surplus

                else:
                    # value of partnerhip for all levels of power
                    Vw_couple = np.zeros(par.num_power)
                    Vm_couple = np.zeros(par.num_power)
                    for iP in range(par.num_power):
                        idx = (t,iP)
                        Vw_couple[iP] = linear_interp.interp_2d(par.grid_love,par.grid_A,sol.Vw_remain_couple[idx],love,A_lag)  # interpolate over love and assets
                        Vm_couple[iP] = linear_interp.interp_2d(par.grid_love,par.grid_A,sol.Vm_remain_couple[idx],love,A_lag)

                    # check participation constraint TODO: should it be the value of REMAINING MARRIED? now it is the value of entering the period as married...
                    Sw = Vw_couple - Vw_single  
                    Sm = Vm_couple - Vm_single
                    power_idx = update_bargaining_index(Sw,Sm,power_idx_lag, par)

                # infer partnership status
                if power_idx < 0.0: # divorce is coded as -1
                    sim.couple[i,t] = False

                else:
                    sim.couple[i,t] = True

            else: # remain single

                sim.couple[i,t] = False

            # update behavior
            if sim.couple[i,t]:
                
                # optimal consumption allocation if couple
                sol_C_tot = sol.C_tot_couple[t,power_idx] 
                C_tot = linear_interp.interp_2d(par.grid_love,par.grid_A,sol_C_tot,love,A_lag)

                sim.Cw_priv[i,t], sim.Cm_priv[i,t], C_pub = intraperiod_allocation(C_tot,power_idx,sol,par)
                sim.Cw_pub[i,t] = C_pub
                sim.Cm_pub[i,t] = C_pub

                # update end-of-period states
                M_resources = resources_couple(A_lag,par) 
                sim.A[i,t] = M_resources - sim.Cw_priv[i,t] - sim.Cm_priv[i,t] - sim.Cw_pub[i,t]
                if t<(par.simT-1):
                    sim.love[i,t+1] = love + par.sigma_love*sim.draw_love[i,t+1]

                # in case of divorce
                sim.Aw[i,t] = par.div_A_share * sim.A[i,t]
                sim.Am[i,t] = (1.0-par.div_A_share) * sim.A[i,t]

                sim.power_idx[i,t] = power_idx
                sim.power[i,t] = par.grid_power[sim.power_idx[i,t]]

            else: # single

                # pick relevant solution for single, depending on whether just became single
                idx_sol_single = t
                sol_single_w = sol.Cw_tot_trans_single[idx_sol_single]
                sol_single_m = sol.Cm_tot_trans_single[idx_sol_single]
                if (power_idx_lag<0):
                    sol_single_w = sol.Cw_tot_single[idx_sol_single]
                    sol_single_m = sol.Cm_tot_single[idx_sol_single]

                # optimal consumption allocations
                Cw_tot = linear_interp.interp_1d(par.grid_Aw,sol_single_w,Aw_lag)
                Cm_tot = linear_interp.interp_1d(par.grid_Am,sol_single_m,Am_lag)
                
                sim.Cw_priv[i,t],sim.Cw_pub[i,t] = intraperiod_allocation_single(Cw_tot,woman,par)
                sim.Cm_priv[i,t],sim.Cm_pub[i,t] = intraperiod_allocation_single(Cm_tot,man,par)

                # update end-of-period states
                Mw = resources_single(Aw_lag,woman,par)
                Mm = resources_single(Am_lag,man,par) 
                sim.Aw[i,t] = Mw - sim.Cw_priv[i,t] - sim.Cw_pub[i,t]
                sim.Am[i,t] = Mm - sim.Cm_priv[i,t] - sim.Cm_pub[i,t]

                # not updated: nans
                # sim.power[i,t] = np.nan
                # sim.love[i,t+1] = np.nan 
                # sim.A[i,t] = np.nan

                sim.power_idx[i,t] = -1

############################
# User-specified functions #
############################
def util(c_priv,c_pub,gender,par,love=0.0):
    if gender == woman:
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

def resources_couple(A,par):
    return par.R*A + par.inc_w + par.inc_m

def resources_single(A,gender,par):
    income = par.inc_w
    if gender == man:
        income = par.inc_m

    return par.R*A + income

def cons_priv_single(C_tot,gender,par):
    # closed form solution for intra-period problem of single.
    if gender == woman:
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

############
# routines #
############
def intraperiod_allocation_single(C_tot,gender,par):
    C_priv = cons_priv_single(C_tot,gender,par)
    C_pub = C_tot - C_priv
    return C_priv,C_pub
 
def intraperiod_allocation(C_tot,iP,sol,par):

    # interpolate pre-computed solution
    j1 = linear_interp.binary_search(0,par.num_Ctot,par.grid_Ctot,C_tot)
    Cw_priv = linear_interp_1d._interp_1d(par.grid_Ctot,sol.pre_Ctot_Cw_priv[iP],C_tot,j1)
    Cm_priv = linear_interp_1d._interp_1d(par.grid_Ctot,sol.pre_Ctot_Cm_priv[iP],C_tot,j1)
    C_pub = C_tot - Cw_priv - Cm_priv 

    return Cw_priv, Cm_priv, C_pub

def solve_intraperiod_couple(C_tot,power,par,starting_val=None):
    
    # setup estimation. Impose constraint that C_tot = Cw+Cm+C
    bounds = optimize.Bounds(0.0, C_tot, keep_feasible=True)
    obj = lambda x: - (power*util(x[0],C_tot-np.sum(x),woman,par) + (1.0-power)*util(x[1],C_tot-np.sum(x),man,par))
    
    # estimate
    x0 = np.array([C_tot/3,C_tot/3]) if starting_val is None else starting_val
    res = optimize.minimize(obj,x0,bounds=bounds)

    # unpack
    Cw_priv = res.x[0]
    Cm_priv = res.x[1]
    C_pub = C_tot - Cw_priv - Cm_priv

    return Cw_priv,Cm_priv,C_pub

def update_bargaining_index(Sw,Sm,iP, par):
    
    # check the participation constraints. Array
    min_Sw = np.min(Sw)
    min_Sm = np.min(Sm)
    max_Sw = np.max(Sw)
    max_Sm = np.max(Sm)

    if (min_Sw >= 0.0) & (min_Sm >= 0.0): # all values are consistent with marriage
        return iP

    elif (max_Sw < 0.0) | (max_Sm < 0.0): # no value is consistent with marriage
        return -1

    else: 
    
        # find lowest (highest) value with positive surplus for women (men)
        Low_w = 0 # in case there is no crossing, this will be the correct value
        Low_m = par.num_power-1 # in case there is no crossing, this will be the correct value
        for _iP in range(par.num_power-1):
            if (Sw[_iP]<0) & (Sw[_iP+1]>=0):
                Low_w = _iP+1
                
            if (Sm[_iP]>=0) & (Sm[_iP+1]<0):
                Low_m = _iP 

        # update the outcomes
        # woman wants to leave
        if iP<Low_w: 
            if Sm[Low_w] > 0: # man happy to shift some bargaining power
                return Low_w
                
            else: # divorce
                return -1
            
        # man wants to leave
        elif iP>Low_m: 
            if Sw[Low_m] > 0: # woman happy to shift some bargaining power
                return Low_m
                
            else: # divorce
                return -1

        else: # no-one wants to leave
            return iP
