import numpy as np
import numba as nb
import scipy.optimize as optimize

from EconModel import EconModelClass
from consav.grids import nonlinspace
from consav import linear_interp, linear_interp_1d
from consav import quadrature
import scipy.stats as stats

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
        self.cpp_filename = '../cppfuncs/solve.cpp'
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
        par.interp_method = 'linear'
        par.precompute_intratemporal = True # if True, precompute intratemporal allocation, else re-solve every time

        par.num_Ctot = 100
        par.max_Ctot = par.max_A*2
        
        par.do_egm = False
        par.num_A_pd = par.num_A * 2
        par.max_A_pd = par.max_A
        par.num_marg_u = 200

        # simulation
        par.seed = 9210
        par.simT = par.T
        par.simN = 50_000
        par.init_A = 0.01
        par.init_love = 0.0
        par.init_power_idx = 10

        # cpp
        par.threads = 8

        par.centered_gradient = True
        
    def allocate(self):
        par = self.par
        sol = self.sol
        sim = self.sim

        # setup grids
        self.setup_grids()
        
        # a. singles
        shape_single = (par.T,par.num_A)                        # single states: T and assets

        # a.1. single to single
        sol.Vw_single_to_single = np.ones(shape_single) - np.inf  # Value
        sol.Vm_single_to_single = np.ones(shape_single) - np.inf

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
        sol.Vw_single_to_single_pd = np.zeros(par.num_A_pd)                # Value of being single, post-decision

        sol.EmargUm_single_to_single_pd = np.zeros(shape_single)          # Expected marginal utility post-decision, man single
        sol.C_totm_single_to_single_pd = np.zeros(par.num_A_pd)           # C for EGM, man single
        sol.Mm_single_to_single_pd = np.zeros(par.num_A_pd)               # Endogenous grid, man single
        sol.Vm_single_to_single_pd = np.zeros(par.num_A_pd)               # Value of being single, post-decision

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

        sol.EVw_cond_meet_partner = np.nan + np.ones(shape_single)
        sol.EVm_cond_meet_partner = np.nan + np.ones(shape_single)


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

        sol.EVw_start_as_couple = np.nan + np.ones(shape_couple)
        sol.EVm_start_as_couple = np.nan + np.ones(shape_couple)
        sol.EmargV_start_as_couple = np.nan + np.ones(shape_couple)
        
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
        sol.pre_Ctot_Cw_priv =  np.nan + np.ones(shape_pre)      # precomputed optimal allocation of consumption over grid of total C 
        sol.pre_Ctot_Cm_priv = np.nan + np.ones(shape_pre)
        sol.pre_Ctot_C_pub = np.nan + np.ones(shape_pre)
    

        # f. simulation
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

        # lifetime utility
        sim.util = np.nan + np.ones((par.simN, par.simT))
        sim.mean_lifetime_util = np.array([np.nan])

        # containers for verifying simulaton
        sim.A_own = np.nan + np.ones(shape_sim)
        sim.A_partner = np.nan + np.ones(shape_sim)

        ## f.1. shocks
        self.allocate_draws()

        ## f.2. initial distribution
        sim.init_A = par.init_A + np.zeros(par.simN)
        sim.init_Aw = sim.init_A * par.div_A_share
        sim.init_Am = sim.init_A * (1.0 - par.div_A_share)
        sim.init_couple = np.ones(par.simN,dtype=np.bool_)
        sim.init_power_idx = par.init_power_idx* np.ones(par.simN,dtype=np.int_)
        sim.init_love = par.init_love + np.zeros(par.simN)
        
        # g. timing
        sol.solution_time = np.array([0.0])

    def allocate_draws(self):
        par = self.par
        sim = self.sim
        shape_sim = (par.simN,par.simT)

        np.random.seed(par.seed)
        sim.draw_love = np.random.normal(size=shape_sim)
        sim.draw_meet = np.random.uniform(size=shape_sim) # for meeting a partner

        sim.draw_uniform_partner_Aw = np.random.uniform(size=shape_sim) # for inverse cdf transformation of partner wealth
        sim.draw_uniform_partner_Am = np.random.uniform(size=shape_sim) # for inverse cdf tranformation of partner wealth

        sim.draw_repartner_love = par.sigma_love*np.random.normal(0.0,1.0,size=shape_sim) #np.random.choice(par.num_love, p=par.prob_partner_love, size=shape_sim) # Love index when repartnering

        
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
        par.grid_marg_u = np.nan + np.ones((par.num_power,par.num_marg_u))
        par.grid_marg_u_for_inv = np.nan + np.ones((par.num_power,par.num_marg_u))

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

        sol = self.sol
        par = self.par 

        # setup grids
        self.setup_grids()

        self.cpp.solve(sol,par)


    def simulate(self):
        sol = self.sol
        sim = self.sim
        par = self.par

        # clear simulation
        for key, val in sim.__dict__.items():
            if 'init' in key or 'draw' in key: continue
            setattr(sim, key, np.zeros(val.shape)+np.nan)

        self.cpp.simulate(sim,sol,par)

        sim.mean_lifetime_util[0] = np.mean(np.sum(sim.util,axis=1))

        

        