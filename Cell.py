# contains class for a single atmospheric layer (Cell), satallites (functions more as a struct), and
# discrete events (Event)

import numpy as np
from BreakupModel import *
from Events import *
import os
import csv
G = 6.67430e-11 # gravitational constant (N*m^2/kg^2)
Me = 5.97219e24 # mass of Earth (kg)
Re = 6371 # radius of Earth (km)

class Cell:
    
    def __init__(self, S_i, D_i, R_i, N_i, logL_edges, chi_edges, event_list, alt, dh, tau_N, v, m_sat, sigma_sat,
                 del_t, lam_sat, alpha_S, alpha_D, alpha_N, alpha_R, P, AM_sat, tau_sat, C_sat, expl_rate_L, 
                 expl_rate_D, m_rb, sigma_rb, lam_rb, AM_rb, tau_rb, C_rb, expl_rate_R):
        '''Constructor for Cell class
    
        Parameter(s):
        S_i : list of initial live satellite values for each satellite type
        D_i : list of initial derelict satellite values for each satellite type
        R_i : list of initial number of rocket bodies of each type
        N_i : initial array of number of debris by L and A/M
        logL_edges : bin edges in log10 of characteristic length (log10(m))
        chi_edges : bin edges in log10(A/M) (log10(m^2/kg))
        event_list : list of discrete events that occur in the cell
        alt : altitude of the shell centre (km)
        dh : width of the shell (km)
        tau_N : 2-d array of atmospheric drag lifetimes for debris (yr)
        v : relative collision speed (km/s)
        m_sat : mass of each satellite type (kg)
        sigma_sat : collision cross-section of each satellite type (m^2)
        del_t : mean satellite lifetime of each type (yr)
        lam_sat : launch rate of satellites of each type into the band (1/yr)
        alpha_S : the fraction of collisions a live satellite of each type fails to avoid with a live satellite
        alpha_D : the fraction of collisions a live satellite of each type fails to avoid with a derelict satellite
        alpha_N : the fraction of collisions a live satellite of each type fails to avoid with trackable debris
        alpha_R : the fraction of collisions a live satellite of each type fails to avoid with a rocket
        P : post-mission disposal probability for satellites of each type
        AM_sat : area-to-mass ratio for satellites of each type (m^2/kg)
        tau_sat : atmospheric drag lifetime of each satellite type (yr)
        C_sat : fit constant for explosions for each satellite type
        expl_rate_L : number of explosions that occur in a 1yr period with a population 
                      of 100 live satellites, for each satellite type
        expl_rate_D : number of explosions that occur in a 1yr period with a population 
                      of 100 derelict satellites, for each satellite type
        m_rb : mass of each rocket body type (kg)
        sigma_rb : collision cross-section of each rocket body type (m^2)
        lam_rb : launch rate of each rocket body type into the shell (1/yr)
        AM_rb : area-to-mass ratio of each rocket body type (m^2/kg)
        tau_rb : atmospheric drag lifetime of each rocket body type (yr)
        C_rb : fit constant for explosions of each rocket body type
        expl_rate_R : number of explosions that occur in a 1yr period with a population 
                      of 100 rocket bodies, for each type

        Output(s):
        Cell instance
        '''

        # setup initial values for tracking satallites
        self.num_sat_types = len(S_i)
        self.num_rb_types = len(R_i)
        self.S = [S_i]
        self.D = [D_i]
        self.m_sat = m_sat
        self.sigma_sat = sigma_sat
        self.sigma_sat_km = self.sigma_sat/1e6 # same thing, but in km^2
        self.del_t = del_t
        self.lam_sat = lam_sat
        self.alpha_S = alpha_S
        self.alpha_D = alpha_D
        self.alpha_R = alpha_R
        self.alpha_N = alpha_N
        self.P = P
        self.AM_sat = AM_sat
        self.tau_sat = tau_sat
        self.C_sat = C_sat
        self.expl_rate_L = expl_rate_L
        self.expl_rate_D = expl_rate_D

        # setup initial values for tracking rockets
        self.R = [R_i]
        self.m_rb = m_rb
        self.sigma_rb = sigma_rb
        self.sigma_rb_km = self.sigma_rb/1e6
        self.lam_rb = lam_rb
        self.AM_rb = AM_rb
        self.tau_rb = tau_rb
        self.C_rb = C_rb
        self.expl_rate_R = expl_rate_R

        # setup initial debris values
        self.N_bins = [N_i]

        # setup other variables
        self.C_c = [0] # catastrophic collisions
        self.C_nc = [0] # non-catastrophic collisions
        self.event_list = event_list
        self.alt = alt
        self.dh = dh
        self.V = 4*np.pi*((6371 + self.alt)**2)*self.dh # volume of the shell
        self.tau_N = tau_N
        self.v = v
        self.v_kyr = self.v*365.25*24*60*60 # convert to km/yr
        self.v_orbit = np.sqrt(G*Me/((Re + alt)*1000))/1000 # orbital velocity in km/s
        self.logL_edges = logL_edges
        self.num_L = len(logL_edges) - 1
        self.logL_ave = np.zeros(self.num_L) # average logL value in each bin
        for i in range(self.num_L):
            self.logL_ave[i] = (logL_edges[i]+logL_edges[i+1])/2
        self.L_ave = 10**self.logL_ave
        self.chi_edges = chi_edges
        self.num_chi = len(chi_edges) - 1
        self.chi_ave = np.zeros(self.num_chi) # average chi value in each bin
        for i in range(self.num_chi):
            self.chi_ave[i] = (chi_edges[i]+chi_edges[i+1])/2
        self.AM_ave = 10**self.chi_ave
        self.trackable = (self.L_ave >= 1/10) # which debris types can be tracked
        self.cat_sat_N = np.full((self.num_sat_types, self.num_L, self.num_chi), False) # tracks which collisions are catastrophic
        self.cat_rb_N = np.full((self.num_rb_types, self.num_L, self.num_chi), False)
        self.update_cat_N()

        # parameters for dxdt calculations
        sigma1_2d = np.resize(self.sigma_sat_km, (self.num_sat_types, self.num_sat_types))
        sigma2_2d = sigma1_2d.transpose()
        self.sigma_comb_satsat = sigma1_2d + sigma2_2d + 2*np.sqrt(sigma1_2d*sigma2_2d) # account for increased cross-section
        self.alphaS1 = np.resize(self.alpha_S, (self.num_sat_types, self.num_sat_types))
        self.alphaS2 = self.alphaS1.transpose()
        self.alphaD1 = np.resize(self.alpha_D, (self.num_sat_types, self.num_sat_types))
        sigma1_2d = np.resize(self.sigma_sat_km, (self.num_rb_types, self.num_sat_types)).transpose()
        sigma2_2d = np.resize(self.sigma_rb_km, (self.num_sat_types, self.num_rb_types))
        self.sigma_comb_satrb = sigma1_2d + sigma2_2d + 2*np.sqrt(sigma1_2d*sigma2_2d) # account for increased cross-section
        self.alphaR1 = np.resize(self.alpha_R, (self.num_rb_types, self.num_sat_types)).transpose()
        sigma1_2d = np.resize(self.sigma_rb_km, (self.num_rb_types, self.num_rb_types))
        sigma2_2d = sigma1_2d.transpose()
        self.sigma_comb_rbrb = sigma1_2d + sigma2_2d + 2*np.sqrt(sigma1_2d*sigma2_2d) # account for increased cross-section
        self.double_count_filter_sat = np.full((self.num_sat_types, self.num_sat_types), False)
        self.double_count_filter_rb = np.full((self.num_rb_types, self.num_rb_types), False)
        for i in range(self.num_sat_types):
            for j in range(i+1,self.num_sat_types):
                self.double_count_filter_sat[i,j] = True
        for i in range(self.num_rb_types):
            for j in range(i+1,self.num_rb_types):
                self.double_count_filter_rb[i,j] = True

    def save(self, filepath, filter, filter_len):
        '''
        saves the current Cell object to .csv and .npz files

        Input(s):
        filepath : explicit path to folder that the files will be saved in (string)
        filter : array of which data points to keep or skip (array of booleans)
        filter_len : number of Trues in the filter

        Output(s): None

        Note(s): event_list is lost, filter should be the same size as the t array from
                 NCell.
        '''

        # save parameters
        csv_file = open(filepath + 'params.csv', 'w', newline='')
        csv_writer = csv.writer(csv_file, dialect='unix')
        csv_writer.writerow([self.num_sat_types, self.num_rb_types, self.alt, self.dh, self.v, self.v_orbit, self.num_L,
                             self.num_chi])
        csv_file.close()

        # write easy arrays
        Cc_array, Cnc_array = np.array(self.C_c)[filter], np.array(self.C_nc)[filter]
        to_save = {'C_c' : Cc_array, 'C_nc' : Cnc_array, 'tau_N' : self.tau_N, 'trackable' : self.trackable,
                   'logL' : self.logL_edges, 'chi' : self.chi_edges}
        np.savez_compressed(filepath + "data.npz", **to_save)

        # write N_bins values
        N_dict = dict()
        index = 0
        for i in range(len(self.N_bins)):
            if filter[i]:
                N_dict[str(index)] = self.N_bins[i]
                index += 1
        np.savez_compressed(filepath + "N_bins.npz", **N_dict)

        # write cat table values
        cat_dict = {'sat' : self.cat_sat_N, 'rb' : self.cat_rb_N}
        np.savez_compressed(filepath + "cat_tables.npz", **cat_dict)

        # write satellites and rockets
        for i in range(self.num_sat_types):
            sat_path = filepath + 'Satellite' + str(i) + '/'
            os.mkdir(sat_path)
            self.save_sat(sat_path, filter, filter_len, i)
        for i in range(self.num_rb_types):
            rb_path = filepath + 'RocketBody' + str(i) + '/'
            os.mkdir(rb_path)
            self.save_rb(rb_path, filter, filter_len, i)

    def save_sat(self, filepath, filter, filter_len, i):
        '''
        saves the current satellite information to .csv and .npz files

        Input(s):
        filepath : explicit path to folder that the files will be saved in (string)
        filter : array of which data points to keep or skip (array of booleans)
        filter_len : number of Trues in the filter
        i : satellite type number

        Output(s): None
        '''

        # save parameters
        csv_file = open(filepath + 'params.csv', 'w', newline='')
        csv_writer = csv.writer(csv_file, dialect='unix')
        csv_writer.writerow([self.m_sat[i], self.sigma_sat[i], self.del_t[i], self.lam_sat[i], 
                             self.alpha_S[i], self.alpha_D[i], self.alpha_N[i], self.alpha_R[i], self.P[i], 
                             self.AM_sat[i], self.tau_sat[i], self.C_sat[i], self.expl_rate_L[i], self.expl_rate_D[i]])
        csv_file.close()

        # save data
        S_array = np.empty(filter_len, dtype=np.double)
        D_array = np.empty(filter_len, dtype=np.double)
        index = 0
        for j in range(len(self.S)):
            if filter[j]:
                S_array[index] = self.S[j][i]
                D_array[index] = self.D[j][i]
                index += 1
        to_save = {'S' : S_array, 'D' : D_array}
        np.savez_compressed(filepath + "data.npz", **to_save)
    
    def save_rb(self, filepath, filter, filter_len, i):
        '''
        saves the current rocket body information to .csv and .npz files

        Input(s):
        filepath : explicit path to folder that the files will be saved in (string)
        filter : array of which data points to keep or skip (array of booleans)
        filter_len : number of Trues in the filter
        i : rocket body type number

        Output(s): None
        '''

        # save parameters
        csv_file = open(filepath + 'params.csv', 'w', newline='')
        csv_writer = csv.writer(csv_file, dialect='unix')
        csv_writer.writerow([self.m_rb[i], self.sigma_rb[i], self.lam_rb[i], self.AM_rb[i], self.tau_rb[i],
                             self.C_rb[i], self.expl_rate_R[i]])
        csv_file.close()

        # save data
        R_array = np.empty(filter_len, dtype=np.double)
        index = 0
        for j in range(len(self.R)):
            if filter[j]:
                R_array[index] = self.R[j][i]
                index += 1
        to_save = {'R' : R_array}
        np.savez_compressed(filepath + "data.npz", **to_save)

    def load(filepath):
        '''
        builds a Cell object from saved data

        Input(s):
        filepath : explicit path to folder that the files are saved in (string)

        Keyword Input(s): None

        Output(s):
        cell : Cell object build from loaded data

        Note(s): cell will not have events
        '''

        cell = Cell.__new__(Cell) # create blank Cell

        # load parameters
        csv_file = open(filepath + 'params.csv', 'r', newline='')
        csv_reader = csv.reader(csv_file, dialect='unix')
        for row in csv_reader: # there's only one row, this extracts it
            cell.num_sat_types = int(row[0])
            cell.num_rb_types = int(row[1])
            cell.alt = float(row[2])
            cell.dh = float(row[3])
            cell.v = float(row[4])
            cell.v_orbit = float(row[5])
            cell.num_L = int(row[6])
            cell.num_chi = int(row[7])
        csv_file.close()

        # load basic arrays
        array_dict = np.load(filepath + "data.npz")
        cell.C_c = array_dict['C_c'].tolist()
        cell.C_nc = array_dict['C_nc'].tolist()
        cell.tau_N = array_dict['tau_N']
        cell.trackable = array_dict['trackable']
        cell.logL_edges = array_dict['logL']
        cell.chi_edges = array_dict['chi']

        # calculate related parameters
        cell.num_L = len(cell.logL_edges) - 1
        cell.logL_ave = np.zeros(cell.num_L) # average logL value in each bin
        for i in range(cell.num_L):
            cell.logL_ave[i] = (cell.logL_edges[i]+cell.logL_edges[i+1])/2
        cell.L_ave = 10**cell.logL_ave
        cell.num_chi = len(cell.chi_edges) - 1
        cell.chi_ave = np.zeros(cell.num_chi) # average logL value in each bin
        for i in range(cell.num_chi):
            cell.chi_ave[i] = (cell.chi_edges[i]+cell.chi_edges[i+1])/2
        cell.AM_ave = 10**cell.chi_ave
        cell.v_kyr = cell.v*365.25*24*60*60 # convert to km/yr
        cell.V = 4*np.pi*((6371 + cell.alt)**2)*cell.dh # volume of the shell

        # load N_bins values
        cell.N_bins = []
        bins_dict = np.load(filepath + "N_bins.npz")
        i = 0
        while True:
            try:
                N_bins = bins_dict[str(i)]
                cell.N_bins.append(N_bins)
            except KeyError:
                break
            i += 1
        
        # load cat table values
        cat_dict = np.load(filepath + "cat_tables.npz")
        cell.cat_sat_N = cat_dict['sat']
        cell.cat_rb_N = cat_dict['rb']

        # setup variables for satellites
        cell.S = []
        cell.D = []
        tot_num_data = len(cell.N_bins) # number of time data points
        for i in range(tot_num_data):
            cell.S.append(np.empty(cell.num_sat_types, dtype=np.double))
            cell.D.append(np.empty(cell.num_sat_types, dtype=np.double))
        cell.m_sat = np.empty(cell.num_sat_types, dtype=np.double)
        cell.sigma_sat = np.empty(cell.num_sat_types, dtype=np.double)
        cell.del_t = np.empty(cell.num_sat_types, dtype=np.double)
        cell.lam_sat = np.empty(cell.num_sat_types, dtype=np.double)
        cell.alpha_S = np.empty(cell.num_sat_types, dtype=np.double)
        cell.alpha_D = np.empty(cell.num_sat_types, dtype=np.double)
        cell.alpha_N = np.empty(cell.num_sat_types, dtype=np.double)
        cell.alpha_R = np.empty(cell.num_sat_types, dtype=np.double)
        cell.P = np.empty(cell.num_sat_types, dtype=np.double)
        cell.AM_sat = np.empty(cell.num_sat_types, dtype=np.double)
        cell.tau_sat = np.empty(cell.num_sat_types, dtype=np.double)
        cell.C_sat = np.empty(cell.num_sat_types, dtype=np.double)
        cell.expl_rate_L = np.empty(cell.num_sat_types, dtype=np.double)
        cell.expl_rate_D = np.empty(cell.num_sat_types, dtype=np.double)

        for i in range(cell.num_sat_types): # load in satellites
            sat_path = filepath + 'Satellite' + str(i) + '/'
            cell.load_sat(sat_path, i)
        
        # compute related parameters
        cell.sigma_sat_km = cell.sigma_sat/1e6

        # setup variables for rockets
        cell.R = []
        for i in range(tot_num_data):
            cell.R.append(np.empty(cell.num_rb_types, dtype=np.double))
        cell.m_rb = np.empty(cell.num_rb_types, dtype=np.double)
        cell.sigma_rb = np.empty(cell.num_rb_types, dtype=np.double)
        cell.lam_rb = np.empty(cell.num_rb_types, dtype=np.double) 
        cell.AM_rb = np.empty(cell.num_rb_types, dtype=np.double)
        cell.tau_rb = np.empty(cell.num_rb_types, dtype=np.double)
        cell.C_rb = np.empty(cell.num_rb_types, dtype=np.double)
        cell.expl_rate_R = np.empty(cell.num_rb_types, dtype=np.double)

        for i in range(cell.num_rb_types):
            rb_path = filepath + 'RocketBody' + str(i) + '/'
            cell.load_rb(rb_path, i)

        # compute related parameters
        cell.sigma_rb_km = cell.sigma_rb/1e6

        cell.event_list = []

        # parameters for dxdt calculations
        sigma1_2d = np.resize(cell.sigma_sat_km, (cell.num_sat_types, cell.num_sat_types))
        sigma2_2d = sigma1_2d.transpose()
        cell.sigma_comb_satsat = sigma1_2d + sigma2_2d + 2*np.sqrt(sigma1_2d*sigma2_2d) # account for increased cross-section
        cell.alphaS1 = np.resize(cell.alpha_S, (cell.num_sat_types, cell.num_sat_types))
        cell.alphaS2 = cell.alphaS1.transpose()
        cell.alphaD1 = np.resize(cell.alpha_D, (cell.num_sat_types, cell.num_sat_types))
        sigma1_2d = np.resize(cell.sigma_sat_km, (cell.num_rb_types, cell.num_sat_types)).transpose()
        sigma2_2d = np.resize(cell.sigma_rb_km, (cell.num_sat_types, cell.num_rb_types))
        cell.sigma_comb_satrb = sigma1_2d + sigma2_2d + 2*np.sqrt(sigma1_2d*sigma2_2d) # account for increased cross-section
        cell.alphaR1 = np.resize(cell.alpha_R, (cell.num_rb_types, cell.num_sat_types)).transpose()
        sigma1_2d = np.resize(cell.sigma_rb_km, (cell.num_rb_types, cell.num_rb_types))
        sigma2_2d = sigma1_2d.transpose()
        cell.sigma_comb_rbrb = sigma1_2d + sigma2_2d + 2*np.sqrt(sigma1_2d*sigma2_2d) # account for increased cross-section
        cell.double_count_filter_sat = np.full((cell.num_sat_types, cell.num_sat_types), False)
        cell.double_count_filter_rb = np.full((cell.num_rb_types, cell.num_rb_types), False)
        for i in range(cell.num_sat_types):
            for j in range(i+1,cell.num_sat_types):
                cell.double_count_filter_sat[i,j] = True
        for i in range(cell.num_rb_types):
            for j in range(i+1,cell.num_rb_types):
                cell.double_count_filter_rb[i,j] = True

        return cell

    def load_sat(self, filepath, i):
        '''
        loads saved satellite information into the current cell

        Input(s):
        filepath : explicit path to folder that the files are saved in (string)
        i : satellite type number

        Output(s): None

        Note(s) : i is a assumed to be a valid number, and that variables are properly initialized
        '''

        # load parameters
        csv_file = open(filepath + 'params.csv', 'r', newline='')
        csv_reader = csv.reader(csv_file, dialect='unix')
        for row in csv_reader: # there's only one row, but this extracts it
            self.m_sat[i], self.sigma_sat[i], self.lam_sat[i] = float(row[0]), float(row[1]), float(row[2])
            self.alpha_S[i], self.alpha_D[i], self.alpha_N[i] = float(row[3]), float(row[4]), float(row[5])
            self.alpha_R[i], self.P[i], self.AM_sat[i] = float(row[6]), float(row[7]), float(row[8])
            self.tau_sat[i], self.C_sat[i], self.expl_rate_L[i] = float(row[9]), float(row[10]), float(row[11])
            self.expl_rate_D[i] = float(row[12])
        csv_file.close()

        # load data
        data_dict = np.load(filepath + "data.npz")
        S_sat = data_dict['S']
        D_sat = data_dict['D']
        for j in range(S_sat.size):
            self.S[j][i] = S_sat[j]
            self.D[j][i] = D_sat[j]

    def load_rb(self, filepath, i):
        '''
        loads saved rocket body information into the current cell

        Input(s):
        filepath : explicit path to folder that the files are saved in (string)
        i : rocket body type number

        Output(s): None

        Note(s) : i is a assumed to be a valid number, and that variables are properly initialized
        '''

        # load parameters
        csv_file = open(filepath + 'params.csv', 'r', newline='')
        csv_reader = csv.reader(csv_file, dialect='unix')
        for row in csv_reader: # there's only one row, but this extracts it
            self.m_rb[i], self.sigma_rb[i], self.lam_rb[i] = float(row[0]), float(row[1]), float(row[2])
            self.AM_rb[i], self.tau_rb[i], self.C_rb[i] = float(row[3]), float(row[4]), float(row[5])
            self.expl_rate_R[i] = float(row[6])
        csv_file.close()

        # load data
        data_dict = np.load(filepath + "data.npz")
        R_rb = data_dict['R']
        for j in range(R_rb.size):
            self.R[j][i] = R_rb[j]

    def dxdt_cell(self, time):
        '''
        calculates the rate of collisions and decays from each debris bin, the rate
        of decaying/de-orbiting satellites, the rate of launches/deorbit starts of satallites, 
        and the rate of creation of derelicts at the given time, due only to events in the cell

        Parameter(s):
        time : index of the values to use

        Keyword Parameter(s): None

        Output(s):
        dSdt : array of rate of change of the number of live satellites in the cell of each type due to only processes
               withing the cell (yr^(-1))
        dDdt : array of rate of change of the number of derelict satellites in the cell of each type
               (excluding derelicts decaying) (yr^(-1))
        dRdt : array of rate of change of number of rocket bodies in the cell of each type (excluding rockets decaying) (yr^(-1))
        D_out : array of rate of satellites decaying from the cell of each type (yr^(-1))
        R_out : array of rate of rocket bodies decaying from the cell of each type (yr^(-1))
        N_out : matrix with the rate of exiting debris from each bin (yr^(-1))
        D_dt : matrix with total rate of collisions between satellites (yr^(-1))
        RD_dt : matrix with total rate of collisions between satellites and rocket bodies (yr^(-1))
        R_dt : matrix with total rate of collisions between rocket bodies (yr^(-1))
        CS_dt : array of matrices with the rate of collisions from each bin with each satellite type (yr^(-1))
        CR_dt : array of matrices with the rate of collisions from each bin with each rocket body type (yr^(-1))
        expl_S : array of rate of explosions for satellites of each type (yr^(-1))
        expl_R : array of rate of explosions for rocket bodies of each type (yr^(-1))

        Note: Assumes that collisions with debris of L_cm < 10cm cannot be avoided, and
        that the given time input is valid
        '''
        
        N_loc = self.N_bins[time]
        S_loc = self.S[time]
        D_loc = self.D[time]
        R_loc = self.R[time]

        # handle satellite-debris collisions
        dSdt = np.resize(N_loc, (self.num_sat_types, self.num_L, self.num_chi)) # collisions with live satallites
        dSdt = np.swapaxes(dSdt, 2, 0) # must be done temporairly for broadcasting
        dDdt = np.array(dSdt) # collisions with derelict satellites
        dSdt *= self.sigma_sat_km*self.v_kyr*S_loc/self.V # compute rates of collision
        dDdt *= self.sigma_sat_km*self.v_kyr*D_loc/self.V
        dSdt[:,self.trackable,:] *= self.alpha_N # account for collision avoidance
        dSdt = np.swapaxes(dSdt,0,2) # switch axes back
        dDdt = np.swapaxes(dDdt,0,2) # switch axes back

        # handle satellite-satellite collisions
        S1 = np.resize(S_loc, (self.num_sat_types, self.num_sat_types))
        S2 = S1.transpose()
        D1 = np.resize(D_loc, (self.num_sat_types, self.num_sat_types))
        D2 = D1.transpose()

        # calculate collision rates
        dSSdt = self.alphaS1*self.alphaS2*self.sigma_comb_satsat*self.v_kyr*S1*S2/self.V
        dSDdt = self.alphaD1*self.sigma_comb_satsat*self.v_kyr*S1*D2/self.V
        dDDdt = self.sigma_comb_satsat*self.v_kyr*D1*D2/self.V  # collisions cannot be avoided

        # compute collisions between satellites and rocket bodies
        S1 = np.resize(S_loc, (self.num_rb_types, self.num_sat_types)).transpose()
        D1 = np.resize(D_loc, (self.num_rb_types, self.num_sat_types)).transpose()
        R2 = np.resize(R_loc, (self.num_sat_types, self.num_rb_types))

        # calculate collision rates
        dSRdt = self.alphaR1*self.sigma_comb_satrb*self.v_kyr*S1*R2/self.V
        dDRdt = self.sigma_comb_satrb*self.v_kyr*D1*R2/self.V # collisions cannot be avoided

        # compute explosion rates for satellites
        expl_S = self.expl_rate_L*S_loc/100
        expl_D = self.expl_rate_D*D_loc/100

        # compute decay/ascend events for satellites
        kill_S, decay_D = S_loc/self.del_t, D_loc/self.tau_sat

        # diagonals are to account for objects of the same type colliding
        dSdt_tot = self.lam_sat - kill_S - np.sum(dSdt, axis=(1,2)) - np.sum(dSSdt, axis=1) - np.sum(dSDdt, axis=1) - np.diag(dSSdt) - np.sum(dSRdt, axis=1) - expl_S
        dDdt_tot = (1-self.P)*kill_S - np.sum(dDdt, axis=(1,2), where=self.cat_sat_N) + np.sum(dSdt, axis=(1,2), where=self.cat_sat_N==False) - np.sum(dSDdt, axis=0) - np.sum(dDDdt, axis=1) - np.diag(dDDdt) - np.sum(dDRdt, axis=1) - expl_D
        CS_dt = dSdt + dDdt # total collisions between satellites and debris

        # handle rocket-debris collisions
        dRdt = np.resize(N_loc, (self.num_rb_types, self.num_L, self.num_chi)) # collisions with rockets
        dRdt = np.swapaxes(dRdt, 0, 2) # must be done temporairly for broadcasting
        dRdt *= self.sigma_rb_km*self.v_kyr*R_loc/self.V # compute rates of collision
        dRdt = np.swapaxes(dRdt,0,2) # switch axes back

        # handle rocket-rocket collisions
        R1 = np.resize(R_loc, (self.num_rb_types, self.num_rb_types))
        R2 = R1.transpose()

        # calculate collision rate
        dRRdt = self.sigma_comb_rbrb*self.v_kyr*R1*R2/self.V

        # handle rocket explosions, decays
        expl_R = self.expl_rate_R*R_loc/100
        decay_R = R_loc/self.tau_rb

        # sum everything up
        dRdt_tot = self.lam_rb - np.sum(dRdt, axis=(1,2), where=self.cat_rb_N) - np.sum(dRRdt, axis=1) - np.diag(dRRdt) - np.sum(dSRdt, axis=0) - np.sum(dDRdt, axis=0) - expl_R

        # calculate decay rates for debris
        decay_N = N_loc/self.tau_N

        # set values to zero to avoid double-counting later on
        dSSdt[self.double_count_filter_sat], dDDdt[self.double_count_filter_sat], dRRdt[self.double_count_filter_rb] = 0, 0, 0

        # compute return values
        D_dt = dSSdt + dSDdt + dDDdt
        RD_dt = dSRdt + dDRdt
        expl_S_tot = expl_S + expl_D

        # return everything
        return dSdt_tot, dDdt_tot, dRdt_tot, decay_D, decay_R, decay_N, D_dt, RD_dt, dRRdt, CS_dt, dRdt, expl_S_tot, expl_R

    def update_cat_N(self):
        '''
        updates values in cat_N based on current mass, v, and bins

        Parameter(s): None

        Keyword Parameter(s): None

        Ouput(s): None
        '''

        for i in range(self.num_sat_types):
            for j in range(self.num_L):
                for k in range(self.num_chi):
                    self.cat_sat_N[i,j,k] = is_catastrophic(self.m_sat[i], self.L_ave[j], self.AM_ave[k], self.v)
        for i in range(self.num_rb_types):
            for j in range(self.num_L):
                for k in range(self.num_chi):
                    self.cat_rb_N[i,j,k] = is_catastrophic(self.m_rb[i], self.L_ave[j], self.AM_ave[k], self.v)
