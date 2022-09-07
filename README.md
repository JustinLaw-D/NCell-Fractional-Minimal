# Minimal Fractional Model

Documentation on the physics model itself is at https://www.overleaf.com/read/wpdhkpfmzxcr.

## General Usage Method

The atmospheric system is represented by an instance of the NCell class. The general method for running a simulation is:

1) Setup an instance of the NCell class to match the desired system, using the \_\_init\_\_ function.
2) Run the simulation of the system using methods of the NCell class.
3) Save the system using the NCell.save() method.
4) Load the system into another Python instance using the NCell.load() static method, and use the getter functions to access data.

## NCell Class
### Initialization
Instances of the NCell class should only be created in two ways, either with the \_\_init\_\_ function or loaded
from a saved object using the NCell.load function. The \_\_init\_\_ function has many parameters, which can be
grouped into the follow categories.

Required Parameters:

S : list of initial number of live satellites in each shell of each type (this is a list of lists, where the first
    index is the cell and the second index is the satellite type)

S_d : list of initial number of deorbiting satellites in each shell of each type (this is a list of lists, where the first
    index is the cell and the second index is the satellite type)

D : list of initial number of derelict satellites in each shell of each type (this is a list of lists, where the first
    index is the cell and the second index is the satellite type)

N_c : initial number of catestrophically lethal debris in each shell (this is a list, where the index is the cell)

alt_edges : edges of the altitude bands to be used (array, km)

lam : launch rate of satellites of each type into each cell (lists of arrays, where the first
      index is the cell and the second index is the satellite type)

The order of indexing for cells is from the cell at the lowest altitude to highest altitude, and the order of satellite
types needs to be the same in all lists. Generally, satellites are considered to be of different types if any of their
parameters (for example mass) differ.

Keyword parameters, which can be broken down into the following categories.

Drag Lifetime Parameters:

update_period : how often the drag lifetimes are updated (yr)

min_lifetime : minimum decay lifetime to allow for debris (yr, this cannot be 0)

CD : drag coefficient for objects

m0 : starting month of the solar cycle (which has 144 months)

min_dt : minimum timestep for calculating decay lifetimes (yr)

max_dt : maximum timestep for calculating decay lifetimes (None or yr)

dtfactor : fraction of altitude/rate of change to take as dt for decay lifetime calculation (yr)

t_max : maximum time to search to for decay lifetime calculation (yr)

setF107 : if not None, value taken for solar flux regardless of current time (None or 10^(-22)W/m^2, default None)

Satellite Parameters:

All of the lists of lists have the first demension as the cell number, and the second the satellite type.
The satellite types need to be in the same order each type, and the cells are always in order of ascending altitude.
All of the 1-D lists are indexed by satellite type.

del_t : mean satellite lifetime for each type (list, yr)

expl_rate_L : number of explosions that occur in a 1yr period with a population of 100 live satellites for
              each type of satellite (list of floats)

expl_rate_D : number of explosions that occur in a 1yr period with a population of 100 derelict satellites
              for each type of satellite (list of floats)

C_sat : fit constant for explosions of each type of satellite (list of floats)

sigma_sat : satellite cross-section of each type (list, m^2, default 10m^2)

alphaS : fraction of collisions with another live satellite that a live satellites of each type fails to 
         avoid in each shell (list of lists)

alphaD : fraction of collisions with another derelict that a live satellites of each type fails to 
         avoid in each shell (list of lists)

alphaN : fraction of collisions with trackable debris that a live satellites of each type fails to 
         avoid in each shell (list of lists)

alphaR : fraction of collisions with a rocket body that a live satellites of each type fails to 
            avoid in each shell (list of lists)

P : post-mission disposal probability for satellites of each type (list)

m_s : mass of the satallites of each type (list, kg)

AM_sat : area-to-mass ratio of the satallites of each type (list, m^2/kg)

Rocket Body Parameters:

All of the lists of lists have the first demension as the cell number, and the second the rocket body type.
The rocket body types need to be in the same order each type, and the cells are always in order of ascending altitude.
All of the 1-D lists are indexed by rocket type.

R_i : list of rocket bodies in each shell of each type (list of lists)

lam_rb : launch rate of rocket bodies of each type into the each shell (list of arrays, 1/yr)

expl_rate_R : number of explosions that occur in a 1yr period with a population of 100 rocket bodies for
              each type of rocket body (list of floats)

C_rb : fit constant for explosions of each type of rocket body (list of floats)

sigma_rb : rocket cross-section of each type (list, m^2)

m_rb : mass of the rocket bodies of each type (list, kg)

AM_rb : area-to-mass ratio of the rocket bodies of each type (list, m^2/kg)

Debris Binning Parameters:

The default values of L_min, L_max, chi_min, and chi_max were tuned to be as restrictive as possible
without losing precision, it's not reccomended that you change them.

L_min : minimum characteristic length to consider (m)

L_max : maximum characteristic length to consider (m)

num_L : number of debris bins in characteristic length

chi_min : minimum log10(A/M) to consider (log10(m^2/kg))

chi_max : maximum log10(A/M) to consider (log10(m^2/kg))

num_chi : number of debris bins in log10(A/M)

delta : initial ratio of the density of disabling to catestrophic debris in each shell

Other Parameters:

The default value of num_dir was tuned to be as restrictive as possible
without losing precision, it's not recommended that you change it. v is calculated
based on the altitude of the shell by default.

num_dir : number of random directions to sample in creating probability tables

table_path : path to saved probability tables (string or None, must be saved in format used in NCell.save)

v : relative collision speed in each shell (list, km/s)

events : the discrete events occuring in the system (list of Event objects, default no events)

### Saving/Loading

The NCell object can be saved/loaded using the corresponding methods. The save method parameters are

filepath : explicit (or relative) path to the directory that the files will be saved in. this directory must exist (string)

name : name of the object, must be a valid unix directory name (string)

gap : largest acceptable time gap between saved data points (yr)

force : whether or not to automatically replace any saved data with the same name (default True)

The adherence to the gap parameter is approximate, so the data may not be saved in the exact time steps you want. The load
function (which is static) is much simpler, having only one parameter

filepath : explicit (or relative) path to the directory that the files are saved in (string)

### Changing Parameters Manually

In general, it's not reccomended that you change any parameters manually. The program is optomized to pre-compute
many values off of them, and so changing them may not have the effect you want. In general, the things it's safe to
change are:

1) Initial values of numbers of satellites, rockets and debris (before running a simulation)
2) Drag lifetime calculation parameters

### Running Simulations

There's two simulation methods built-in to the NCell class. The first is an Euler method (run_sim_euler), which has parameters

T : length of the simulation (yr)

dt : timestep used by the simulation (yr, it's reccomended that you set this to 1/1000 at most)

upper : whether or not to have debris come into the top shell (bool, default True)

The second is a predictor-corrector method with adaptive timestep (run_sim_precor) which has parameters

T : length of the simulation (yr)

dt_i : initial timestep used by the simulation (yr)

dt_min : minimum time step used by the simulation is (yr)

dt_max : maximum time step used by simulation (yr)

tolerance : tolerance for adaptive time step

err_factor : how close to tolerance epsilon can be without actually triggering a redo

upper : whether or not to have debris come into the top shell (bool, default True)

A few notes about both simulations. First, if your system includes altitudes below 500km, it's reccomended that you set
dt_min (or dt in the case of an Euler simulation) to at most the minimum decay lifetime, as decay lifetimes tend to hit
the limit (which is 1e-4yr by default). Second, if your system starts off being relatively empty, it's reccomended that you
set the initial time step to a small value (as the tolerance checker igonres parameters with a value of zero). Third, the adaptive time step in the predictor-corrector method makes no attempt to check how accurate the number of collisions estimation is, and won't adjust the time step based on it.

### Getter Methods

The NCell class also contains 8 getter methods. They can be grouped into the following categories.

1) Time getter (get_t). Gets the list of time (in yr) that data points are saved at.
2) Satellite/Rocket getters (get_S, get_D, get_R). These get the number of live and de-orbiting satellites,
and the number of rocket bodies respectively. They are indexed by (cell, type, time), with the cells in ascending order.
3) Debris getter (get_N). This gets the total amount of debris in each cell at each time, and is indexed by (cell, time), with cells in ascending order.
4) Collision getters (get_C, get_Cc, get_Cnc). These get the total, catestrophic, and non-catestrophic collisions respectively in each cell at each time, indexed in the same way as get_N.

## Event Class

The Event class is essentially a template class, designed to be overriden. It represents discrete events (such as a anti-satellite weapons test) that could occur in the system. In general, it has three parameters:

alt : altitude of the event (km)

time : list of times that the event occurs (yr or None)

freq : frequency the event occurs at (1/yr or None)

at least one of time or freq must not be None. The class also has a function, run_event, with the following signature

Input(s):

S : number of live satellites of each type in the current cell (list of floats)

D : number of derelict satellites of each type in the current cell (list of floats)

R : number of rocket bodies of each type in the current cell (list of floats)

N : binned amount of debris in current cell (2d array)

logL_edges : logL edge values for the debris bins (log10(m))

chi_edges : chi edge values for the debris bins (log10(m^2/kg))

Output(s):

dS : change in the number of live satellites of each type in the current cell (array of floats)

dD : change in the number of derelict satellites of each type in the cell (array of floats)

dR : change in the number of rocket bodies of each type in the cell (array of floats)

dN : change in the number of debris in the curren cell, not including debris
     produced by collisions

coll : list of collisions occuring in the current cell in the form [(kg, kg, typ, #)],
       i.e. [(m1, m2, typ, number of collisions)]. typ can be one of 'sat' (satellite-satellite),
       'sr' (satellite-rocket, where satellite is m1), or 'rb' (rocket-rocket)

expl : list of explosions occuring in the current cell in the form [(C, typ, #)], where
       C is the relevant fit constant and typ is the type of body exploding ('sat' or 'rb)

By default the return values are all zero and empty list. In creating an event, you should create your own event class which inherits from event and customizes the run_event function (although the run_event function must still have the same signature). Two such classes are provided in Events.py, and are

1) ExplEvent, which represents a group of explosions occuring, and takes the additional parameter expl_list (the list of explosions that run_event will return) in it's \_\_init\_\_ function.
2) CollEvent, which is the same thing for collisions.
