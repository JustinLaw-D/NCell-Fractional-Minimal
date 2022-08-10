# contains class for collection of cells representing orbital shells

from Cell import *
from Events import *
import numpy as np
from BreakupModel import *
from AtmosphericDecayModels import *
from copy import deepcopy
import os
import shutil
import csv

G = 6.67430e-11 # gravitational constant (N*m^2/kg^2)
Re = 6371 # radius of Earth (km)
Me = 5.97219e24 # mass of Earth (kg)

class NCell:

    def __init__(self, S, D, N_c, alt_edges, lam, update_period=1/12, min_lifetime=0, CD=2.2, m0=0, min_dt=0, 
                 max_dt=0.1, dtfactor=1/100, t_max=np.inf, setF107=None, events=[], R_i=None, lam_rb=None, del_t=None, 
                 expl_rate_L=None, expl_rate_D=None, C_sat=None, sigma_sat=None, expl_rate_R=None, C_rb=None, sigma_rb=None, 
                 v=None, delta=None, alphaS=None, alphaD=None, alphaN=None, alphaR=None, P=None, m_s=None, m_rb=None, AM_sat=None, 
                 AM_rb=None, L_min=1e-3, L_max=1, num_L=10, chi_min=-2, chi_max=1.0, num_chi=10, num_dir=1000, table_path=None):
        '''
        Constructor for NCell class
    
        Parameter(s):
        S : list of initial number of live satellites in each shell of each type (list of arrays)
        D : list of initial number of derelict satellites in each shell of each type (list of arrays)
        N_c : initial number of catestrophically lethal debris in each shell (array)
        alt_edges : edges of the altitude bands to be used (array, km)
        lam : launch rate of satellites of each type into each shell (list of arrays, 1/yr)

        Keyword Parameter(s):
        update_period : how often the drag lifetimes are updated (yr, default 1/12)
        min_lifetime : minimum decay lifetime to allow for debris (yr, default 0)
        CD : drag coefficient for objects (default 2.2)
        m0 : starting month of the solar cycle (default 0)
        min_dt : minimum timestep for calculating decay lifetimes (yr, default 0)
        max_dt : maximum timestep for calculating decay lifetimes (None or yr, default 0.1)
        dtfactor : fraction of altitude/rate of change to take as dt for decay lifetime calculation (yr, default 1/100)
        t_max : maximum time to search to for decay lifetime calculation (yr, default infinite)
        setF107 : if not None, value taken for solar flux regardless of current time (None or 10^(-22)W/m^2, default None)
        events : the discrete events occuring in the system (list of Event objects, default no events)
        R_i : list of rocket bodies in each shell of each type (list of lists, default no rocket bodies)
        lam_rb : launch rate of rocket bodies of each type into the each shell (list of arrays, 1/yr, default all 0)
        del_t : mean satellite lifetime of each type (list, yr, default 5yr)
        expl_rate_L : number of explosions that occur in a 1yr period with a population of 100 live satellites for
                      each type of satellite (list of floats, default all 0)
        expl_rate_D : number of explosions that occur in a 1yr period with a population of 100 derelict satellites
                      for each type of satellite (list of floats, default expl_rate_L)
        C_sat : fit constant for explosions of each type of satellite (list of floats, default all 1)
        sigma_sat : satellite cross-section of each type (list, m^2, default 10m^2)
        expl_rate_R : number of explosions that occur in a 1yr period with a population of 100 rocket bodies for
                      each type of rocket body (list of floats, default all 0)
        C_rb : fit constant for explosions of each type of rocket body (list of floats, default all 1)
        sigma_rb : rocket cross-section of each type (list, m^2, default 10m^2)
        v : relative collision speed in each shell (list, km/s, default 10km/s)
        delta : initial ratio of the density of disabling to catestrophic debris in each shell (list, default 10)
        alphaS : fraction of collisions with another live satellite that a live satellites of each type fails to 
                 avoid in each shell (list of lists, default 0)
        alphaD : fraction of collisions with another derelict that a live satellites of each type fails to 
                 avoid in each shell (list of lists, default alphaN)
        alphaN : fraction of collisions with trackable debris that a live satellites of each type fails to 
                 avoid in each shell (list of lists, default 0.2)
        alphaR : fraction of collisions with a rocket body that a live satellites of each type fails to 
                 avoid in each shell (list of lists, default alphaD)
        P : post-mission disposal probability for satellites of each type (list, default 0.95)
        m_s : mass of the satallites of each type (list, kg, default 250kg)
        m_rb : mass of the rocket bodies of each type (list, kg, default 250kg)
        AM_sat : area-to-mass ratio of the satallites of each type (list, m^2/kg, default 1/(20*2.2)m^2/kg)
        AM_rb : area-to-mass ratio of the rocket bodies of each type (list, m^2/kg, default 1/(20*2.2)m^2/kg)
        L_min : minimum characteristic length to consider (m, default 1mm)
        L_max : maximum characteristic length to consider (m, default 1m)
        num_L : number of debris bins in characteristic length (default 10)
        chi_min : minimum log10(A/M) to consider (log10(m^2/kg), default -2)
        chi_max : maximum log10(A/M) to consider (log10(m^2/kg), default 1)
        num_chi : number of debris bins in log10(A/M) (default 10)
        num_dir : number of random directions to sample in creating probability tables (default 1000)
        table_path : path to saved probability tables (string or None, must be saved in format used in NCell.save)

        Output(s):
        NCell instance

        Note: no size checks are done on the arrays, the program will crash if any of the arrays are of incorrect size.
        shells are assumed to be given in order of ascending altitude. if you only want to pass values in the
        keyword argument for certain shells, put None in the list for all other shells. internally, cells may have
        padded space in their arrays, use the getter functions to clean those up.
        '''

        if len(S) == 0: # check if there's no shells
            print("ERROR: System must have at least one shell!")
            return
        self.num_cells = len(S) # total number of cells in the system

        # convert Nones to array of Nones
        if events is None:
            events = []
        if R_i is None:
            R_i = [[]]*self.num_cells
        if lam_rb is None:
            lam_rb = [None]*self.num_cells
        if v is None:
            v = [10]*self.num_cells
        if delta is None:
            delta = [10]*self.num_cells
        if alphaS is None:
            alphaS = [None]*self.num_cells
        if alphaD is None:
            alphaD = [None]*self.num_cells
        if alphaN is None:
            alphaN = [None]*self.num_cells
        if alphaR is None:
            alphaR = [None]*self.num_cells

        self.num_sat_types = len(S[0]) # save number of satellite and rocket types
        self.num_rb_types = len(R_i[0])
        self.update_period = update_period # save lifetime calculation parameters
        self.min_lifetime = min_lifetime
        self.CD = CD
        self.m0 = m0
        self.min_dt = min_dt
        self.max_dt = max_dt
        self.dtfactor = dtfactor
        self.t_max = t_max
        self.setF107 = setF107
        self.alts = np.zeros(len(alt_edges)-1) # setup altitude bins
        self.dhs = np.zeros(self.alts.shape)
        for i in range(self.alts.size):
            self.dhs[i] = alt_edges[i+1]-alt_edges[i]
            self.alts[i] = (alt_edges[i]+alt_edges[i+1])/2
        self.num_L = num_L
        self.num_chi = num_chi
        self.time = 0 # index of current time step
        self.lupdate_time = 0 # index of last time drag lifetimes were updated
        self.t = [0] # list of times traversed
        self.cells = [] # start list of cells
        # generate bins for log10(L), chi
        self.logL_edges = np.linspace(np.log10(L_min), np.log10(L_max), num=num_L+1)
        self.logL_ave = np.zeros(self.num_L) # average logL value in each bin
        for i in range(self.num_L):
            self.logL_ave[i] = (self.logL_edges[i]+self.logL_edges[i+1])/2
        self.L_ave = 10**self.logL_ave
        self.chi_edges = np.linspace(chi_min, chi_max, num=num_chi+1)
        self.chi_ave = np.zeros(self.num_chi) # average chi value in each bin
        for i in range(self.num_chi):
            self.chi_ave[i] = (self.chi_edges[i]+self.chi_edges[i+1])/2
        self.AM_ave = 10**self.chi_ave
        self.bin_masses = np.empty((self.num_L, self.num_chi)) # average mass in each bin
        for i in range(self.num_L):
            A = find_A(self.L_ave[i])
            for j in range(self.num_chi):
                self.bin_masses[i,j] = A/self.AM_ave[j]
        self.num_dir = num_dir

        for i in range(self.num_cells): # iterate through shells

            # convert Nones to array of Nones
            if lam_rb[i] is None:
                lam_rb[i] = [None]*self.num_rb_types
            if del_t is None:
                del_t = [None]*self.num_sat_types
            if expl_rate_L is None:
                expl_rate_L = [None]*self.num_sat_types
            if expl_rate_D is None:
                expl_rate_D = [None]*self.num_sat_types
            if C_sat is None:
                C_sat = [None]*self.num_sat_types
            if sigma_sat is None:
                sigma_sat = [None]*self.num_sat_types
            if expl_rate_R is None:
                expl_rate_R = [None]*self.num_rb_types
            if C_rb is None:
                C_rb = [None]*self.num_rb_types
            if sigma_rb is None:
                sigma_rb = [None]*self.num_rb_types
            if alphaS[i] is None:
                alphaS[i] = [None]*self.num_sat_types
            if alphaD[i] is None:
                alphaD[i] = [None]*self.num_sat_types
            if alphaN[i] is None:
                alphaN[i] = [None]*self.num_sat_types
            if alphaR[i] is None:
                alphaR[i] = [None]*self.num_sat_types
            if P is None:
                P = [None]*self.num_sat_types
            if m_s is None:
                m_s = [None]*self.num_sat_types
            if m_rb is None:
                m_rb = [None]*self.num_rb_types
            if AM_sat is None:
                AM_sat = [None]*self.num_sat_types
            if AM_rb is None:
                AM_rb = [None]*self.num_rb_types

            S_cell = np.empty(self.num_sat_types, dtype=np.double) # setup satellite parameter arrays
            D_cell = np.empty(self.num_sat_types, dtype=np.double)
            m_sat_cell = np.empty(self.num_sat_types, dtype=np.double)
            sigma_sat_cell = np.empty(self.num_sat_types, dtype=np.double)
            del_t_cell = np.empty(self.num_sat_types, dtype=np.double)
            lam_sat_cell = np.empty(self.num_sat_types, dtype=np.double)
            alpha_S_cell = np.empty(self.num_sat_types, dtype=np.double)
            alpha_D_cell = np.empty(self.num_sat_types, dtype=np.double)
            alpha_R_cell = np.empty(self.num_sat_types, dtype=np.double)
            alpha_N_cell = np.empty(self.num_sat_types, dtype=np.double)
            P_cell = np.empty(self.num_sat_types, dtype=np.double)
            AM_sat_cell = np.empty(self.num_sat_types, dtype=np.double)
            tau_sat_cell = np.empty(self.num_sat_types, dtype=np.double)
            C_sat_cell = np.empty(self.num_sat_types, dtype=np.double)
            expl_rate_L_cell = np.empty(self.num_sat_types, dtype=np.double)
            expl_rate_D_cell = np.empty(self.num_sat_types, dtype=np.double)

            for j in range(self.num_sat_types): # iterate through satellite types, and generate parameters for each
                
                # convert Nones to default values
                if del_t[j] is None:
                    del_t[j] = 5
                if expl_rate_L[j] is None:
                    expl_rate_L[j] = 0
                if expl_rate_D[j] is None:
                    expl_rate_D[j] = expl_rate_L[j]
                if C_sat[j] is None:
                    C_sat[j] = 1
                if sigma_sat[j] is None:
                    sigma_sat[j] = 10
                if alphaS[i][j] is None:
                    alphaS[i][j] = 0
                if alphaN[i][j] is None:
                    alphaN[i][j] = 0.2
                if alphaD[i][j] is None:
                    alphaD[i][j] = alphaN[i][j]
                if alphaR[i][j] is None:
                    alphaR[i][j] = alphaD[i][j]
                if P[j] is None:
                    P[j] = 0.95
                if m_s[j] is None:
                    m_s[j] = 250
                if AM_sat[j] is None:
                    AM_sat[j] = 1/(20*2.2)

                # compute atmospheric drag lifetime for satallites in the shell
                tau = max(drag_lifetime(self.alts[i] + self.dhs[i]/2, self.alts[i] - self.dhs[i]/2, AM_sat[j], CD, 1/365.25, m0,
                                        min_dt, max_dt, dtfactor, t_max, setF107), self.min_lifetime)
                S_cell[j] = S[i][j]
                D_cell[j] = D[i][j]
                m_sat_cell[j] = m_s[j]
                sigma_sat_cell[j] = sigma_sat[j]
                del_t_cell[j] = del_t[j]
                lam_sat_cell[j] = lam[i][j]
                alpha_S_cell[j] = alphaS[i][j]
                alpha_D_cell[j] = alphaD[i][j]
                alpha_R_cell[j] = alphaR[i][j]
                alpha_N_cell[j] = alphaN[i][j]
                P_cell[j] = P[j]
                AM_sat_cell[j] = AM_sat[j]
                tau_sat_cell[j] = tau
                C_sat_cell[j] = C_sat[j]
                expl_rate_L_cell[j] = expl_rate_L[j]
                expl_rate_D_cell[j] = expl_rate_D[j]

            # setup rocket-body parameter arrays
            R_cell = np.empty(self.num_rb_types, dtype=np.double) # setup satellite parameter arrays
            lam_rb_cell = np.empty(self.num_rb_types, dtype=np.double)
            m_rb_cell = np.empty(self.num_rb_types, dtype=np.double)
            sigma_rb_cell = np.empty(self.num_rb_types, dtype=np.double)
            AM_rb_cell = np.empty(self.num_rb_types, dtype=np.double)
            tau_rb_cell = np.empty(self.num_rb_types, dtype=np.double)
            C_rb_cell = np.empty(self.num_rb_types, dtype=np.double)
            expl_rate_R_cell = np.empty(self.num_rb_types, dtype=np.double)

            for j in range(self.num_rb_types): # iterate through rocket types, and generate object for each
                
                # convert Nones to default values
                if lam_rb[i][j] is None:
                    lam_rb[i][j] = 0
                if expl_rate_R[j] is None:
                    expl_rate_R[j] = 0
                if C_rb[j] is None:
                    C_rb[j] = 1
                if sigma_rb[j] is None:
                    sigma_rb[j] = 10
                if m_rb[j] is None:
                    m_rb[j] = 250
                if AM_rb[j] is None:
                    AM_rb[j] = 1/(20*2.2)

                # compute atmospheric drag lifetime for rocket bodies in the shell
                tau = max(drag_lifetime(self.alts[i] + self.dhs[i]/2, self.alts[i] - self.dhs[i]/2, AM_rb[j], CD, 1/365.25, m0,
                                        min_dt, max_dt, dtfactor, t_max, setF107), self.min_lifetime)
                R_cell[j] = R_i[i][j]
                lam_rb_cell[j] = lam_rb[i][j]
                m_rb_cell[j] = m_rb[j]
                sigma_rb_cell[j] = sigma_rb[j]
                AM_rb_cell[j] = AM_rb[j]
                tau_rb_cell[j] = tau
                C_rb_cell[j] = C_rb[j]
                expl_rate_R_cell[j] = expl_rate_R[j]

            # calculate decay paremeters for debris, initial debris values
            N_initial, tau_N = np.zeros((num_L, num_chi)), np.empty(num_chi, dtype=np.double)
            # generate initial distributions
            for j in range(self.num_L):
                bin_L = 0
                bin_bot_L, bin_top_L = self.logL_edges[j], self.logL_edges[j+1]
                if (bin_bot_L < -1) and (bin_top_L > -1):
                    lam_factor = (-1-bin_bot_L)/(bin_top_L-bin_bot_L)
                    bin_L += lam_factor*N_c[i]*delta[i]*(L_cdf(1e-1, L_min, 1e-1, 'expl') - L_cdf(10**bin_bot_L, L_min, 1e-1, 'expl'))
                    bin_L += (1-lam_factor)*N_c[i]*(L_cdf(10**bin_top_L, 1e-1, L_max, 'expl') - L_cdf(1e-1, 1e-1, L_max, 'expl'))
                elif bin_bot_L >= -1:
                    bin_L += N_c[i]*(L_cdf(10**bin_top_L, 1e-1, L_max, 'expl') - L_cdf(10**bin_bot_L, 1e-1, L_max, 'expl'))
                else:
                    bin_L += N_c[i]*delta[i]*(L_cdf(10**bin_top_L, L_min, 1e-1, 'expl') - L_cdf(10**bin_bot_L, L_min, 1e-1, 'expl'))
                N_initial[j,0] = bin_L # put everything in the lowest A/M bin
            for j in range(self.num_chi):
                tau_N[j] = max(drag_lifetime(self.alts[i] + self.dhs[i]/2, self.alts[i] - self.dhs[i]/2, self.AM_ave[j], CD, 1/365.25, 
                                             m0, min_dt, max_dt, dtfactor, t_max, setF107), self.min_lifetime)

            # figure out which events are in this cell
            events_loc = []
            for event in events:
                if (event.alt > self.alts[i] - self.dhs[i]/2) and (event.alt <= self.alts[i] + self.dhs[i]/2) : events_loc.append(event)

            # initialize cell
            cell = Cell(S_cell, D_cell, R_cell, N_initial, self.logL_edges, self.chi_edges, events_loc,
                        self.alts[i], self.dhs[i], tau_N, v[i], m_sat_cell, sigma_sat_cell, del_t_cell, lam_sat_cell, 
                        alpha_S_cell, alpha_D_cell, alpha_N_cell, alpha_R_cell, P_cell, AM_sat_cell, 
                        tau_sat_cell, C_sat_cell, expl_rate_L_cell, expl_rate_D_cell, m_rb_cell, sigma_rb_cell, lam_rb_cell, 
                        AM_rb_cell, tau_rb_cell, C_rb_cell, expl_rate_R_cell)
            self.cells.append(cell)
            if i == self.num_cells - 1: self.upper_N = deepcopy(N_initial) # take the debris field above to be initial debris of top

        # generate uniformly distributed directions using Fibbonacci spiral
        phi, theta = np.empty(num_dir), np.empty(num_dir)
        golden = (1+np.sqrt(5))/2 # golden ratio
        for i in range(num_dir):
            x = (i/golden) % 1
            y = i/num_dir
            phi[i] = 2*np.pi*x
            theta[i] = np.arccos(1-2*y)
        # check how probability tables will be aquired
        if table_path is None: # generate tables
            # setup probability tables (arguments are collision bin, final bin, logL, chi)
            self.sat_coll_probability_tables = np.zeros((self.num_cells, self.num_cells, self.num_L, self.num_chi))
            self.rb_coll_probability_tables = np.zeros((self.num_cells, self.num_cells, self.num_L, self.num_chi))
            self.sat_expl_probability_tables = np.zeros((self.num_cells, self.num_cells, self.num_L, self.num_chi))
            self.rb_expl_probability_tables = np.zeros((self.num_cells, self.num_cells, self.num_L, self.num_chi))
            self.fill_prob_tables(phi, theta) # compute probability tables
        else: # load tables
            prob_dict = np.load(table_path)
            self.sat_coll_probability_tables = prob_dict['sat_coll_tables']
            self.rb_coll_probability_tables = prob_dict['rb_coll_tables']
            self.sat_expl_probability_tables = prob_dict['sat_expl_tables']
            self.rb_expl_probability_tables = prob_dict['rb_expl_tables']

    def fill_prob_tables(self, phi, theta):
        '''
        calculates probability tables

        Input(s):
        phi : array of phi components of directions
        theta : array of theta components of directions

        Keyword Input(s): None

        Output(s): None
        '''

        L_min, L_max = 10**self.logL_edges[0], 10**self.logL_edges[-1] # edges of the parameter space
        chi_min, chi_max = self.chi_edges[0], self.chi_edges[-1]
        for i in range(self.num_cells): # iterate through where the event occurs
            v0 = self.cells[i].v_orbit*1000 # orbital velocity in m/s
            r = self.cells[i].alt # in km
            for j in range(self.num_cells): # iterate through final location cells
                curr_cell = self.cells[j]
                alt_min = curr_cell.alt - curr_cell.dh/2 # in km
                alt_max = curr_cell.alt + curr_cell.dh/2
                v_min2 = G*Me*(2/((Re + r)*1000) - 1/((Re + alt_min)*1000)) # minimum velocity squared (m/s)
                v_max2 = G*Me*(2/((Re + r)*1000) - 1/((Re + alt_max)*1000)) # maximum velocity squared (m/s)
                # handle vprime_cdf
                if v_min2 < 0 and v_max2 < 0 : pass
                else:
                    if v_min2 < 0:
                        coll_probs = vprime_cdf(np.sqrt(v_max2), v0, theta, phi, self.chi_ave, 'coll')
                        expl_probs = vprime_cdf(np.sqrt(v_max2), v0, theta, phi, self.chi_ave, 'expl')
                    else:
                        coll_probs = vprime_cdf(np.sqrt(v_max2), v0, theta, phi, self.chi_ave, 'coll') - vprime_cdf(np.sqrt(v_min2), v0, theta, phi, self.chi_ave, 'coll')
                        expl_probs = vprime_cdf(np.sqrt(v_max2), v0, theta, phi, self.chi_ave, 'expl') - vprime_cdf(np.sqrt(v_min2), v0, theta, phi, self.chi_ave, 'expl')
                    # do monte-carlo integration
                    sum_coll = np.sum(coll_probs, 1)
                    sum_expl = np.sum(expl_probs, 1)
                    self.sat_coll_probability_tables[i,j,:,:] = sum_coll/self.num_dir # save the results
                    self.rb_coll_probability_tables[i,j,:,:] = sum_coll/self.num_dir
                    self.sat_expl_probability_tables[i,j,:,:] = sum_expl/self.num_dir
                    self.rb_expl_probability_tables[i,j,:,:] = sum_expl/self.num_dir

        # probability of L being in each bin
        L_prob_coll = L_cdf(10**self.logL_edges[1:], L_min, L_max, 'coll') - L_cdf(10**self.logL_edges[:-1], L_min, L_max, 'coll')
        L_prob_expl = L_cdf(10**self.logL_edges[1:], L_min, L_max, 'expl') - L_cdf(10**self.logL_edges[:-1], L_min, L_max, 'expl')
        # probability of chi, for each bin
        chi_prob_sat = X_cdf(self.chi_edges[1:], chi_min, chi_max, self.L_ave, 'sat') - X_cdf(self.chi_edges[:-1], chi_min, chi_max, self.L_ave, 'sat')
        chi_prob_rb = X_cdf(self.chi_edges[1:], chi_min, chi_max, self.L_ave, 'rb') - X_cdf(self.chi_edges[:-1], chi_min, chi_max, self.L_ave, 'rb')
        # compute total result
        for i in range(self.num_chi): # this loop is needed or broadcasting gets ugly
            self.sat_coll_probability_tables[:,:,:,i] *= L_prob_coll
            self.rb_coll_probability_tables[:,:,:,i] *= L_prob_coll
            self.sat_expl_probability_tables[:,:,:,i] *= L_prob_expl
            self.rb_expl_probability_tables[:,:,:,i] *= L_prob_expl
        self.sat_coll_probability_tables *= chi_prob_sat
        self.sat_expl_probability_tables *= chi_prob_sat
        self.rb_coll_probability_tables *= chi_prob_rb
        self.rb_expl_probability_tables *= chi_prob_rb

    def save(self, filepath, name, gap=0, force=True):
        '''
        saves the current NCell object to .csv and .npz files

        Input(s):
        filepath : explicit path to folder that the files will be saved in (string)
        name : name of the object, must be a valid unix folder name (string)

        Keyword Input(s):
        gap : largest acceptable time gap between saved data points (yr, default 0 i.e. save all data)
        force : whether or not to automatically replace any saved data with the same name (default True)

        Output(s): None

        Note(s): events are lost. adherence to the "gap" value is approximate, and may behave
        strangely if the time step is close to the gap size.
        '''

        true_path = filepath + name + '/'
        try:
            os.mkdir(true_path) # make the folder representing the object
        except FileExistsError:
            x = 'y'
            if not force:
                x = input("File with this name already exists. Replace it (y/n): ")
            if x == 'y':
                shutil.rmtree(true_path)
                os.mkdir(true_path)
            else : return

        # write parameters
        csv_file = open(true_path + 'params.csv', 'w', newline='')
        csv_writer = csv.writer(csv_file, dialect='unix')
        if self.setF107 == None:
            write_F = -1 # use -1 to mean None
        else:
            write_F = self.setF107
        if self.max_dt == None:
            write_maxdt = -1
        else:
            write_maxdt = self.max_dt
        if self.t_max == np.inf:
            write_tmax = -1
        else:
            write_tmax = self.dt_max
        csv_writer.writerow([self.num_L, self.num_chi, self.num_cells, self.num_dir, self.min_lifetime, self.CD,
                             self.m0, self.min_dt, write_maxdt, self.dtfactor, write_tmax, write_F])
        csv_file.close()

        # write easy arrays
        t_arr = np.array(self.t)
        filter = np.full(t_arr.shape, False) # build filter based on time steps
        filter_len = 0 # number of Trues in the filter
        if t_arr.size > 0:
            prev_t = t_arr[0]
            filter[0] = True
            filter_len += 1
            for i in range(1, t_arr.size):
                if t_arr[i] - prev_t >= gap:
                    prev_t = t_arr[i]
                    filter[i] = True
                    filter_len += 1
        to_save = {'alts' : self.alts, 'dhs' : self.dhs, 't' : t_arr[filter], 'logL' : self.logL_edges, 'chi' : self.chi_edges}
        np.savez_compressed(true_path + "data.npz", **to_save)
        
        # save probability tables
        to_save = {'sat_coll_tables' : self.sat_coll_probability_tables, 'sat_expl_tables' : self.sat_expl_probability_tables,
                   'rb_coll_tables' : self.rb_coll_probability_tables, 'rb_expl_tables' : self.rb_expl_probability_tables}
        np.savez_compressed(true_path + "prob_tables.npz", **to_save)

        # save the Cells
        for i in range(self.num_cells):
            cell_path = true_path + "cell" + str(i) + "/"
            os.mkdir(cell_path)
            self.cells[i].save(cell_path, filter, filter_len)

    def load(filepath):
        '''
        builds an NCell object from saved data

        Input(s):
        filepath : explicit path to folder that the files are saved in (string)

        Keyword Input(s): None

        Output(s):
        atmos : NCell object build from loaded data

        Note(s): atmos will not have events
        '''

        atmos = NCell.__new__(NCell) # empty initialization

        # load parameters
        csv_file = open(filepath + 'params.csv', 'r', newline='')
        csv_reader = csv.reader(csv_file, dialect='unix')
        for row in csv_reader: # there's only one row, this extracts it
            atmos.num_L = int(row[0])
            atmos.num_chi = int(row[1])
            atmos.num_cells = int(row[2])
            atmos.num_dir = int(row[3])
            atmos.min_lifetime = float(row[4])
            atmos.CD = float(row[5])
            atmos.m0 = float(row[6])
            atmos.min_dt = float(row[7])
            atmos.max_dt = float(row[8])
            if atmos.max_dt == -1 : atmos.max_dt = None
            atmos.dtfactor = float(row[9])
            atmos.t_max = float(row[10])
            if atmos.t_max == -1 : atmos.t_max = np.inf
            atmos.setF107 = float(row[11])
            if atmos.setF107 == -1 : atmos.setF107 = None
        csv_file.close()

        # load in simple numpy arrays
        array_dict = np.load(filepath + 'data.npz')
        atmos.alts = array_dict['alts']
        atmos.dhs = array_dict['dhs']
        atmos.t = array_dict['t'].tolist()
        atmos.time = len(atmos.t) - 1 # set time to the end of the data
        atmos.lupdate_time = atmos.time
        atmos.logL_edges = array_dict['logL']
        atmos.chi_edges = array_dict['chi']

        # compute related parameters
        atmos.logL_ave = np.zeros(atmos.num_L) # average logL value in each bin
        for i in range(atmos.num_L):
            atmos.logL_ave[i] = (atmos.logL_edges[i]+atmos.logL_edges[i+1])/2
        atmos.L_ave = 10**atmos.logL_ave
        atmos.chi_ave = np.zeros(atmos.num_chi)
        for i in range(atmos.num_chi):
            atmos.chi_ave[i] = (atmos.chi_edges[i]+atmos.chi_edges[i+1])/2
        atmos.AM_ave = 10**atmos.chi_ave
        atmos.bin_masses = np.empty((atmos.num_L, atmos.num_chi)) # average mass in each bin
        for i in range(atmos.num_L):
            A = find_A(atmos.L_ave[i])
            for j in range(atmos.num_chi):
                atmos.bin_masses[i,j] = A/atmos.AM_ave[j]

        # load in probability tables
        prob_dict = np.load(filepath + "prob_tables.npz")
        atmos.sat_coll_probability_tables = prob_dict['sat_coll_tables']
        atmos.rb_coll_probability_tables = prob_dict['rb_coll_tables']
        atmos.sat_expl_probability_tables = prob_dict['sat_expl_tables']
        atmos.rb_expl_probability_tables = prob_dict['rb_expl_tables']

        # get Cells
        atmos.cells = []
        for i in range(atmos.num_cells):
            cell_path = filepath + "cell" + str(i) + "/"
            atmos.cells.append(Cell.load(cell_path))
        atmos.num_sat_types = len(atmos.cells[0].m_sat)
        atmos.num_rb_types = len(atmos.cells[0].m_rb)

        return atmos

    def dxdt(self, time, upper):
        '''
        calculates the rates of change of all parameters at the given time

        Parameter(s):
        time : time (index) of the values to be used
        upper : whether or not to have debris come into the top shell (bool)

        Keyword Parameter(s): None

        Output(s):
        dSdt : list of rates of change in S for each cell (2-d array, 1/yr)
        dDdt : list of rates of change in D for each cell (2-d array, 1/yr)
        dRdt : list of rates of change in R for each cell (2-d array, 1/yr)
        dNdt : list of rates of change in the N matrix for each cell (3-d array, 1/yr)
        dCcdt : list of rates of change in C_c for each cell (array, 1/yr)
        dCcldt : list of rates of change in C_nc for each cell (array1/yr)

        Note : does not check that the time input is valid. arrays are indexed in order (cell, type)
        '''

        top_cell = self.cells[-1]
        top_Nin = self.upper_N/top_cell.tau_N # debris going into top cell
        dSdt = np.zeros((self.num_cells, self.num_sat_types)) # array of changes in satallite values
        dDdt = np.zeros((self.num_cells, self.num_sat_types)) # array of changes in derelict values
        dRdt = np.zeros((self.num_cells, self.num_rb_types)) # array of changes in rocket body values
        dNdt =  np.zeros((self.num_cells, self.num_L, self.num_chi)) # array of changes in debris values
        sat_coll =  np.zeros((self.num_cells, self.num_sat_types, self.num_sat_types)) # array of satellite-satellite collisions
        RS_coll = np.zeros((self.num_cells, self.num_sat_types, self.num_rb_types)) # array of rocket-satellite collisions
        R_coll = np.zeros((self.num_cells, self.num_rb_types, self.num_rb_types)) # array of rocket-rocket collisions
        NS_coll = np.zeros((self.num_cells, self.num_sat_types, self.num_L, self.num_chi)) # array of collision values for satellites
        NR_coll = np.zeros((self.num_cells, self.num_rb_types, self.num_L, self.num_chi)) # array of collision values for rockets
        NS_expl = np.zeros((self.num_cells, self.num_sat_types)) # array of explosion values for satellites
        NR_expl = np.zeros((self.num_cells, self.num_rb_types)) # array of explosion values for rockets

        # get initial D_in, N_in values
        D_in = np.zeros((self.num_cells+1, self.num_sat_types))
        R_in = np.zeros((self.num_cells+1, self.num_rb_types))
        N_in  = np.zeros((self.num_cells+1, self.num_L, self.num_chi))
        if upper : N_in[-1,:,:] = top_Nin

        # iterate through cells, from top to bottom
        for i in range(self.num_cells):
            curr_cell = self.cells[i]
            dSdt[i,:], dDdt[i,:], dRdt[i,:], D_in[i,:], R_in[i,:], N_in[i,:,:], sat_coll[i,:,:], RS_coll[i,:,:], R_coll[i,:,:], NS_coll[i,:,:,:], NR_coll[i,:,:,:], NS_expl[i,:], NR_expl[i,:] = curr_cell.dxdt_cell(time)
            # simulate collisions and explosions
            self.sim_colls(dNdt, sat_coll[i,:,:], curr_cell.m_sat, curr_cell.m_sat, i, 'sat') # sat-sat
            self.sim_colls_satrb(dNdt, RS_coll[i,:,:], curr_cell.m_sat, i, 'sat') # sat-rb
            self.sim_colls_satrb(dNdt, RS_coll[i,:,:], curr_cell.m_rb, i, 'rb') # sat-rb
            self.sim_colls(dNdt, NS_coll[i,:,:,:], curr_cell.m_sat, self.bin_masses, i, 'sat') # sat-debris
            self.sim_expl(dNdt, NS_expl, curr_cell.C_sat, i, 'sat') # sat explosions
            self.sim_colls(dNdt, R_coll[i,:,:], curr_cell.m_rb, curr_cell.m_rb, i, 'rb') # rb-rb
            self.sim_colls(dNdt, NR_coll[i,:,:,:], curr_cell.m_rb, self.bin_masses, i, 'rb') # rb-debris
            self.sim_expl(dNdt, NR_expl, curr_cell.C_rb, i, 'rb') # rb explosions
                    
        # add on debris lost to collisions
        if self.num_sat_types != 0:
            dNdt -= np.sum(NS_coll, axis=1)
        if self.num_rb_types != 0:
            dNdt -= np.sum(NR_coll, axis=1)

        # go through cells from bottom to top to correct values
        dDdt += D_in[1:,:] - D_in[:self.num_cells,:]
        dRdt += R_in[1:,:] - R_in[:self.num_cells,:]
        dNdt += N_in[1:,:,:] - N_in[:self.num_cells,:,:]

        # update values
        dCcdt = np.sum(sat_coll, axis=(1,2)) + np.sum(RS_coll, axis=(1,2)) + np.sum(R_coll, axis=(1,2))
        dCncdt = np.zeros(self.num_cells)
        for i in range(self.num_cells):
            curr_cell = self.cells[i]
            dCcdt[i] += np.sum(NS_coll[i,:,:,:][curr_cell.cat_sat_N]) + np.sum(NR_coll[i,:,:,:][curr_cell.cat_rb_N])
            dCncdt[i] += np.sum(NS_coll[i,:,:,:][curr_cell.cat_sat_N==False]) + np.sum(NR_coll[i,:,:,:][curr_cell.cat_rb_N==False])

        return dSdt, dDdt, dRdt, dNdt, dCcdt, dCncdt

    def run_sim_euler(self, T, dt=1, upper=True):
        '''
        simulates the evolution of the debris-satallite system for T years using a Euler method

        Parameter(s):
        T : length of the simulation (yr)

        Keyword Parameter(s):
        dt : timestep used by the simulation (yr, default 1yr)
        upper : whether or not to have debris come into the top shell (bool, default True)

        Output(s): None
        '''

        self.sim_events() # run initial discrete events

        while self.t[self.time] < T:
            if (self.t[self.time] - self.t[self.lupdate_time]) >= self.update_period:
                    self.update_lifetimes(self.t[self.time])
                    self.lupdate_time = self.time
            dSdt, dDdt, dRdt, dNdt, dCcdt, dCncdt = self.dxdt(self.time, upper) # get current rates of change

            for i in range(self.num_cells): # iterate through cells and update values
                curr_cell = self.cells[i]
                curr_cell.S.append(curr_cell.S[self.time] + dSdt[i]*dt)
                curr_cell.D.append(curr_cell.D[self.time] + dDdt[i]*dt)
                curr_cell.R.append(curr_cell.R[self.time] + dRdt[i]*dt)
                curr_cell.N_bins.append(curr_cell.N_bins[self.time] + dNdt[i]*dt)
                curr_cell.C_c.append(curr_cell.C_c[self.time] + dCcdt[i]*dt)
                curr_cell.C_nc.append(curr_cell.C_nc[self.time] + dCncdt[i]*dt)
            self.t.append(self.t[self.time] + dt) # update time
            self.time += 1
            self.sim_events() # run discrete events

    def run_sim_precor(self, T, dt_i=1, dt_min=1/1000, dt_max=1, tolerance=1, err_factor=1e-6, upper=True):
        '''
        simulates the evolution of the debris-satallite system for T years using predictor-corrector model

        Parameter(s):
        T : length of the simulation (yr)

        Keyword Parameter(s):
        dt_i : initial timestep used by the simulation (yr, default 1)
        dt_min : minimum time step used by the simulation is (yr, default 1/1000)
        dt_max : maximum time step used by simulation (yr, default 1)
        tolerance : tolerance for adaptive time step
        err_factor : how close to tolerance epsilon can be without actually triggering a redo
        upper : whether or not to have debris come into the top shell (bool, default True)

        Output(s): None

        Note(s): AB(2) method is used as predictor, Trapezoid method as corrector
        '''

        warning_given = False # whether or not a warning has been given yet
        # get additional initial value if needed
        if self.time == 0 : self.run_sim_euler(dt_min, dt=dt_min, upper=upper)
        # get previous rate of change values
        self.update_lifetimes(self.t[self.time-1])
        dSdt_n, dDdt_n, dRdt_n, dNdt_n, dCcdt_n, dCncdt_n = self.dxdt(self.time-1, upper=upper)
        # get current rate of change values
        self.update_lifetimes(self.t[self.time])
        self.lupdate_time = self.time
        dSdt_n1, dDdt_n1, dRdt_n1, dNdt_n1, dCcdt_n1, dCncdt_n1 = self.dxdt(self.time, upper=upper)
        dt_old = self.t[self.time] - self.t[self.time-1] # set up old time step variable
        dt = dt_i
        updated, redo = False, False

        while self.t[self.time] < T:
            if updated and redo:
                self.update_lifetimes(self.t[self.time])
            elif updated:
                self.lupdate_time = self.time
            redo = False
            updated = False
            # step forwards using AB(2) method
            for i in range(self.num_cells): # iterate through cells and update values
                curr_cell = self.cells[i]

                if len(curr_cell.N_bins) < self.time + 2: # check if we need to lengthen things
                    curr_cell.S.append(0)
                    curr_cell.D.append(0)
                    curr_cell.R.append(0)
                    curr_cell.N_bins.append(0)
                    curr_cell.C_c.append(0)
                    curr_cell.C_nc.append(0)

                curr_cell.S[self.time+1] = curr_cell.S[self.time] + 0.5*dt*((2+dt/dt_old)*dSdt_n1[i]-(dt/dt_old)*dSdt_n[i])
                curr_cell.D[self.time+1] = curr_cell.D[self.time] + 0.5*dt*((2+dt/dt_old)*dDdt_n1[i]-(dt/dt_old)*dDdt_n[i])
                curr_cell.R[self.time+1] = curr_cell.R[self.time] + 0.5*dt*((2+dt/dt_old)*dRdt_n1[i]-(dt/dt_old)*dRdt_n[i])
                curr_cell.N_bins[self.time+1] = curr_cell.N_bins[self.time] + 0.5*dt*((2+dt/dt_old)*dNdt_n1[i]-(dt/dt_old)*dNdt_n[i])
                curr_cell.C_c[self.time+1] = curr_cell.C_c[self.time] + 0.5*dt*((2+dt/dt_old)*dCcdt_n1[i]-(dt/dt_old)*dCcdt_n[i])
                curr_cell.C_nc[self.time+1] = curr_cell.C_nc[self.time] + 0.5*dt*((2+dt/dt_old)*dCncdt_n1[i]-(dt/dt_old)*dCncdt_n[i])
            # get predicted rate of change from AB(2) method prediction
            if (self.t[self.time] + dt - self.t[self.lupdate_time]) >= self.update_period:
                self.update_lifetimes(self.t[self.time] + dt)
                updated = True
            dSdt_n2, dDdt_n2, dRdt_n2, dNdt_n2, dCcdt_n2, dCncdt_n2 = self.dxdt(self.time+1, upper=upper)
            # set up variable for step size checking
            epsilon = 0
            # re-do step using Trapezoid method
            for i in range(self.num_cells): # iterate through cells and update values
                curr_cell = self.cells[i]
                # get old values
                old_S = curr_cell.S[self.time+1]
                old_D = curr_cell.D[self.time+1]
                old_R = curr_cell.R[self.time+1]
                old_N = curr_cell.N_bins[self.time+1]

                # update with new values
                curr_cell.S[self.time+1] = curr_cell.S[self.time] + 0.5*(dSdt_n2[i]+dSdt_n1[i])*dt
                curr_cell.D[self.time+1] = curr_cell.D[self.time] + 0.5*(dDdt_n2[i]+dDdt_n1[i])*dt
                curr_cell.R[self.time+1] = curr_cell.R[self.time] + 0.5*(dRdt_n2[i]+dRdt_n1[i])*dt
                curr_cell.N_bins[self.time+1] = curr_cell.N_bins[self.time] + 0.5*(dNdt_n2[i]+dNdt_n1[i])*dt

                # estimate errors with old and new values
                valid_choice_S = curr_cell.S[self.time] != 0
                valid_choice_D = curr_cell.D[self.time] != 0
                valid_choice_R = curr_cell.R[self.time] != 0
                valid_choice_N = curr_cell.N_bins[self.time] != 0
                if np.any(valid_choice_S) == True:
                    epsilon_options = np.abs((1/3)*(dt/(dt+dt_old))*(curr_cell.S[self.time+1][valid_choice_S]-old_S[valid_choice_S]))
                    epsilon = max(np.amax(epsilon_options), epsilon)
                if np.any(valid_choice_D) == True:
                    epsilon_options = np.abs((1/3)*(dt/(dt+dt_old))*(curr_cell.D[self.time+1][valid_choice_D]-old_D[valid_choice_D]))
                    epsilon = max(np.amax(epsilon_options), epsilon)
                if np.any(valid_choice_R) == True:
                    epsilon_options = np.abs((1/3)*(dt/(dt+dt_old))*(curr_cell.R[self.time+1][valid_choice_R]-old_R[valid_choice_R]))
                    epsilon = max(np.amax(epsilon_options), epsilon)
                if np.any(valid_choice_N) == True:
                    epsilon_options = np.abs((1/3)*(dt/(dt+dt_old))*(curr_cell.N_bins[self.time+1][valid_choice_N]-old_N[valid_choice_N]))
                    epsilon = max(np.amax(epsilon_options), epsilon)
   
                # we don't really care that much about the accuracy of the collision count
                curr_cell.C_c[self.time+1] = curr_cell.C_c[self.time] + 0.5*(dCcdt_n2[i]+dCcdt_n1[i])*dt
                curr_cell.C_nc[self.time+1] = curr_cell.C_nc[self.time] + 0.5*(dCncdt_n2[i]+dCncdt_n1[i])*dt

            # update step size, and check if calculation needs to be redone
            if (epsilon > tolerance) and (np.abs(epsilon - tolerance) > err_factor):
                redo = True
            new_dt = min(np.abs(dt*(tolerance/epsilon)**(1/3)), dt_max)
            if redo:
                if dt <= dt_min:
                    if not warning_given:
                        print('WARNING : System may be too stiff to integrate')
                        warning_given = True
                    redo=False
                    new_dt = dt_min
                else:
                    dt = new_dt
                    continue

            # update time
            self.t.append(self.t[self.time] + dt)
            self.time += 1
            dt_old = dt
            dt = new_dt
            # run events
            self.sim_events()
            # update which are the old and new rates of change
            dSdt_n, dDdt_n, dRdt_n, dNdt_n, dCcdt_n, dCncdt_n = dSdt_n1, dDdt_n1, dRdt_n1, dNdt_n1, dCcdt_n1, dCncdt_n1
            dSdt_n1, dDdt_n1, dRdt_n1, dNdt_n1, dCcdt_n1, dCncdt_n1 = self.dxdt(self.time, upper)

    def sim_colls(self, dNdt, rate, m_1, m_2, indx, typ):
        '''
        updates dNdt by distributing rates of collisions between two objects of mass m_1, m_2 in
        the index'th cell
        
        Parameter(s):
        dNdt : current dNdt values (3-d array, 1/yr)
        rate : rate of collisions to simulate (2-d or 3-d array, 1/yr)
        m_1 : mass of the first object (array, kg)
        m_2 : mass of the second object (1-d or 2-d array, kg)
        indx : index of the cell the collision occurs in
        typ : object type of the objects, either 'sat' (satellite) or 'rb' (rocket body)

        Keyword Parameter(s): None

        Output(s): if rate is 2-d, m_2 must be 1-d. if rate is 3-d, m_2 must be 2-d
        '''
        v_rel = self.cells[indx].v # collision velocity (km/s)
        M = calc_M(m_1, m_2, v_rel) # M factor
        Lmin, Lmax = 10**self.logL_edges[0], 10**self.logL_edges[-1] # min and max characteristic lengths
        N_debris = np.sum(calc_Ntot(M, Lmin, Lmax, 'coll')*rate) # total rate of debris creation
        if typ == 'sat':
            dNdt += self.sat_coll_probability_tables[indx,:,:,:]*N_debris
        elif typ == 'rb':
            dNdt += self.rb_coll_probability_tables[indx,:,:,:]*N_debris

    def sim_colls_satrb(self, dNdt, rate, m, indx, typ):
        '''
        version of sim_coll used for the satellite-rocket body collisions workaround, where
        each object is simulated as having its own catastrophic collision
        
        Parameter(s):
        dNdt : current dNdt values (3-d array, 1/yr)
        rate : rates of collisions to simulate (2-d array, 1/yr)
        m : mass of the object (array, kg)
        indx : index of the cell the collision occurs in
        typ : object type in the collision, either 'sat' (satellite) or 'rb' (rocket body)

        Keyword Parameter(s): None

        Output(s): None

        Note(s): rate must be indexed in order (satellite, rocket body)
        '''

        Lmin, Lmax = 10**self.logL_edges[0], 10**self.logL_edges[-1] # min and max characteristic lengths
        if typ == 'sat':
            # total rate of debris creation, we sum over num_rb_types to get debris produced for each type of collision
            N_debris = np.sum(calc_Ntot(m, Lmin, Lmax, 'coll')*np.sum(rate, axis=1))
            dNdt += self.sat_coll_probability_tables[indx,:,:,:]*N_debris
        elif typ == 'rb':
            # total rate of debris creation, we sum over num_sat_types to get debris produced for each type of collision
            N_debris = np.sum(calc_Ntot(m, Lmin, Lmax, 'coll')*np.sum(rate,axis=0))
            dNdt += self.rb_coll_probability_tables[indx,:,:,:]*N_debris

    def sim_expl(self, dNdt, rate, C, indx, typ):
        '''
        updates dNdt by distributing a rate of explosions for an object with constant C in
        the index'th cell
        
        Parameter(s):
        dNdt : current dNdt values (3-d matrix, 1/yr)
        rate : rate of explosions to simulate (array, 1/yr)
        C : fit constant for the explosion (can be array)
        indx : index of the cell the collision occurs in
        typ : object type of the main object, either 'sat' (satellite) or 'rb' (rocket body)

        Keyword Parameter(s): None

        Output(s): None
        '''

        Lmin, Lmax = 10**self.logL_edges[0], 10**self.logL_edges[-1] # min and max characteristic lengths
        N_debris = np.sum(calc_Ntot(0, Lmin, Lmax, 'expl', C=C)*rate) # total rate of debris creation
        if typ == 'sat':
            dNdt += self.sat_expl_probability_tables[indx,:,:,:]*N_debris
        elif typ == 'rb':
            dNdt += self.rb_expl_probability_tables[indx,:,:,:]*N_debris

    def sim_events(self):
        '''
        simulates discrete events at the current time

        Input(s): None

        Keyword Input(s): None

        Output(s): None
        '''

        dN = np.zeros((self.num_cells, self.num_L, self.num_chi)) # debris change matrix

        for i in range(self.num_cells):

            curr_cell = self.cells[i]
            dS, dD = np.zeros(curr_cell.num_sat_types), np.zeros(curr_cell.num_sat_types)
            dR = np.zeros(curr_cell.num_rb_types)
            dN_loc = np.zeros((self.num_L, self.num_chi)) # debris change from non-collision sources
            coll_list = []
            expl_list = []
            S, D, R = curr_cell.S[self.time], curr_cell.D[self.time], curr_cell.R[self.time]
            N = curr_cell.N_bins[self.time]

            for event in curr_cell.event_list: # iterate through possible events

                if event.time is not None: # events at specific times
                    while event.time != [] and (event.time[0] <= self.t[self.time]):
                        dS_temp, dD_temp, dR_temp, dN_loc_temp, coll_temp, expl_temp = event.run_event(S, D, R, N, self.logL_edges, self.chi_edges)
                        event.time.pop(0)
                        dS += dS_temp
                        dD += dD_temp
                        dR += dR_temp
                        dN_loc += dN_loc_temp
                        coll_list.extend(coll_temp)
                        expl_list.extend(expl_temp)

                if event.freq is not None: # events occuring at specific frequencies
                    if self.t[self.time] - event.last_event >= 1/event.freq:
                        dS_temp, dD_temp, dR_temp, dN_loc_temp, coll_temp, expl_temp = event.run_event(S, D, R, N, self.logL_edges, self.chi_edges)
                        dS += dS_temp
                        dD += dD_temp
                        dR += dR_temp
                        dN_loc += dN_loc_temp
                        coll_list.extend(coll_temp)
                        expl_list.extend(expl_temp)
                        event.last_event = self.t[self.time]

                # update values
                curr_cell.S[self.time] += dS
                curr_cell.D[self.time] += dD
                curr_cell.R[self.time] += dR
                curr_cell.N_bins[self.time] += dN_loc

                # handle collisions and explosions
                self.parse_coll(dN, coll_list, i)
                self.parse_expl(dN, expl_list, i)

        # update with debris from collisions/explosions
        for i in range(self.num_cells):
            curr_cell = self.cells[i]
            curr_cell.N_bins[self.time] += dN[i,:,:]

    def parse_coll(self, dN, coll_list, i):
        '''
        parses and runs discrete collision events, storing the debris generated in dN

        Input(s):
        dN : 3d matrix of changes in debris for each bin and cell
        coll_list : list of collisions occuring in the current cell in the form [(kg, kg, typ, #)],
                    i.e. [(m1, m2, typ, number of collisions)]. typ can be one of 'sat' (satellite-satellite),
                    'sr' (satellite-rocket, where satellite is m1), or 'rb' (rocket-rocket)
        i : index of the current cell
        '''

        for coll in coll_list: # iterate through list
            m1, m2, typ, num = coll # unpack the list
            m1, m2, num = np.array(m1), np.array(m2), np.reshape(np.array(num), (1,1)) # needed for compatibility
            if typ == 'sat' or typ == 'rb':
                self.sim_colls(dN, num, m1, m2, i, typ)
            elif typ == 'sr':
                self.sim_colls_satrb(dN, num, m1, i, 'sat')
                self.sim_colls_satrb(dN, num, m2, i, 'rb')
    
    def parse_expl(self, dN, expl_list, i):
        '''
        parses and runs discrete explosion events, storing the debris generated in dN

        Input(s):
        dN : 3d matrix of changes in debris for each bin and cell
        expl_list : list of explosions occuring in the current cell in the form [(C, typ, #)], where
                    C is the relevant fit constant and typ is the type of body exploding ('sat' or 'rb)
        i : index of the current cell
        '''

        for expl in expl_list: # iterate through list
            C, typ, num = expl # unpack the list
            C, num = np.array(C), np.array(num) # needed for compatibility
            self.sim_expl(dN, num, C, i, typ)

    def update_lifetimes(self, t):
        '''
        updates all drag lifetimes in the system

        Input(s):
        t : time to call drag_lifetime at (yr)

        Keyword Input(s): None

        Output(s): None
        '''

        for i in range(self.num_cells): # iterate through cells
            curr_cell = self.cells[i]
            alt = curr_cell.alt
            dh = curr_cell.dh
            for j in range(self.num_sat_types): # handle satellites
                curr_cell.tau_sat[j] = max(drag_lifetime(alt+dh/2, alt-dh/2, curr_cell.AM_sat[j], self.CD, 1/365.25, self.m0 + t*12, 
                                                         self.min_dt, self.max_dt, self.dtfactor, self.t_max, self.setF107), self.min_lifetime)
            for j in range(self.num_rb_types): # handle rockets
                curr_cell.tau_rb[j] = max(drag_lifetime(alt+dh/2, alt-dh/2, curr_cell.AM_rb[j], self.CD, 1/365.25, self.m0 + t*12, 
                                                        self.min_dt, self.max_dt, self.dtfactor, self.t_max, self.setF107), self.min_lifetime)
            for j in range(self.num_chi): # handle debris
                curr_cell.tau_N[j] = max(drag_lifetime(alt+dh/2, alt-dh/2, curr_cell.AM_ave[j], self.CD, 1/365.25, self.m0 + t*12, 
                                                       self.min_dt, self.max_dt, self.dtfactor, self.t_max, self.setF107), self.min_lifetime)

    def get_t(self):
        '''
        returns array of times used in the simulation

        Parameter(s): None

        Keyword Parameter(s): None

        Returns:
        array of t values (yr)
        '''

        return self.t
    
    def get_S(self):
        '''
        returns list of lists of lists for number of live satellites in each shell

        Parameter(s): None

        Keyword Parameter(s): None

        Returns:
        list of array of S values for each cell of each type, in order of ascending altitude
        '''

        to_return = []
        for cell in self.cells:
            to_return.append([])
            for j in range(cell.num_sat_types):
                to_return[-1].append([])
                for k in range(self.time+1):
                    to_return[-1][j].append(cell.S[k][j])
        return to_return

    def get_D(self):
        '''
        returns list of lists of lists for number of derelict satellites in each shell

        Parameter(s): None

        Keyword Parameter(s): None

        Returns:
        list of array of D values for each cell of each type, in order of ascending altitude
        '''

        to_return = []
        for cell in self.cells:
            to_return.append([])
            for j in range(cell.num_sat_types):
                to_return[-1].append([])
                for k in range(self.time+1):
                    to_return[-1][j].append(cell.D[k][j])
        return to_return

    def get_R(self):
        '''
        returns list of lists of lists for number of rocket bodies in each shell

        Parameter(s): None

        Keyword Parameter(s): None

        Returns:
        list of array of R values for each cell of each type, in order of ascending altitude
        '''

        to_return = []
        for cell in self.cells:
            to_return.append([])
            for j in range(cell.num_rb_types):
                to_return[-1].append([])
                for k in range(self.time+1):
                    to_return[-1][j].append(cell.R[k][j])
        return to_return

    def get_N(self):
        '''
        returns arrays for number of debris in each shell

        Parameter(s): None

        Keyword Parameter(s): None

        Returns:
        list of array of total N values for each cell, in order of ascending altitude
        '''

        to_return = []
        for cell in self.cells:
            N = []
            for i in range(len(cell.N_bins)):
                N.append(np.sum(cell.N_bins[i]))
            to_return.append(np.array(N))
        return to_return

    def get_C(self):
        '''
        returns arrays for number of collisions in each shell

        Parameter(s): None

        Keyword Parameter(s): None

        Returns:
        list of array of total C values for each cell, in order of ascending altitude
        '''

        to_return = []
        for cell in self.cells:
            to_return.append(np.array(cell.C_c) + np.array(cell.C_nc))
        return to_return
    
    def get_Cc(self):
        '''
        returns arrays for number of catastrophic collisions in each shell

        Parameter(s): None

        Keyword Parameter(s): None

        Returns:
        list of array of C_c values for each cell, in order of ascending altitude
        '''

        to_return = []
        for cell in self.cells:
            to_return.append(cell.C_c)
        return to_return

    def get_Cnc(self):
        '''
        returns arrays for number of non-catastrophic collisions in each shell

        Parameter(s): None

        Keyword Parameter(s): None

        Returns:
        list of array of C_nc values for each cell, in order of ascending altitude
        '''

        to_return = []
        for cell in self.cells:
            to_return.append(cell.C_nc)
        return to_return

    def alt_to_index(self, h):
        '''
        Converts given altitude to cell index

        Parameter(s):
        h : altitude to convert (km)

        Keyword Parameter(s): None

        Output(s):
        index : index corresponding to that altitude, or -1 if none is found
        '''

        for i in range(self.num_cells):
            alt, dh = self.alts[i], self.dh[i]
            if (alt - dh/2 <= h) and (alt + dh/2 >= h):
                return i
        return -1
