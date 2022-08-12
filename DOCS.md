# DOCS

This document contains more thorough documentation of the finer details of the code. It's reccomended that you read this before attempting to manually change parameters, or if you're having issues with the code.

## Cell Class

Starting with the cell class, we'll break the member variables down into three categories. For specifics on what the variables represent, see the README or the docstring of Cell.\_\_init\_\_().

### Safe to change
The following variables can be safely changed by the user, without much worry of breaking the program. The condition that the variable match the required shape (i.e. right shape for the number of satellite types, rocket body types, and debris bins) is always assumed.

1) alpha_S/alpha_D/alpha_N/alpha_R
2) lam_sat/lam_rb

### Safe to change (with percausions)
The following variables can be changed manually by the user, as long as certain conditions are met.

1) S/S_d/R/N_bins : All of these can have their initial values changed (i.e the first values in the list). The requirements are that all values must be positive, N_bins must match the given binning shape, the other initial values must match the number of satellite and rocket body types respectively, and the program cannot have run a simulation yet. Changing any later values in the list would probably be fine under the same restrictions, as long as you didn't try simulating with the predictor-corrector method afterwords, as it uses a range of past data points for future predictions. I don't reccomend it.
2) event_list : This can safely be changed at any time, as long as the lupdate_time parameter of the Event (see events section of this documentation) is properly set up.
3) del_t/expl_rate_L/expl_rate_D/expl_rate_R/C_sat/C_rb : All of these can be safely changed but must agree across all cells.
4) sigma_sat/sigma_rb : These can be changed, but under three conditions. First, their values must agree across all cells. Second, the values of sigma_sat_km/sigma_rb_km (the same variables but in km^2 instead of m^2) must be updated in each cell. Third, the arrays sigma_comb_satsat, sigma_comb_satrb, and sigma_comb_rbrb must all be updated in each cell. See the code of Cell.\_\_init\_\_() for how to update (i.e. create) these arrays.
5) alpha_S/alpha_D/alpha_N/alpha_R : All of these can be changed safely, as long as the arrays alphaS1, alphaS2, alphaD1, and alphaR1 are updated. See the code of Cell.\_\_init\_\_() for how to update (i.e. create) these arrays.
6) v : This can be safely changed, as long as v_kyr (the same value in km/yr instead of km/s) is updated as well, and the method Cell.update_cat_N() is called after updating the values.
7) m_s/m_rb : These values can be safely changed under two conditions. First, their values agree across all cells. Second, the method Cell.update_cat_N() must be called after updating the values
8) AM_sat/AM_rb : These values can be safely changed under two conditions. First, their values agree across all cells. Second, the after updating the values, run the following bit of code, assuming your NCell object is called system
```Python3
system.update_lifetimes(system.t[system.time])
system.lupdate_time = system.time
```

### Unsafe to change
Any variables not mentioned above could, in theory, be safely changed by the user, but this can only be garunteed if you're familiar with the entirety of the code.

## NCell Class

Again for the NCell class, we'll break the member variables down into three categories. For specifics on what the variables represent, see the README or the docstring of NCell.\_\_init\_\_().

### Safe to change
The following variables can be safely changed by the user, without much worry of breaking the program. The condition that the variable match the required shape (i.e. right shape for the number of satellite types, rocket body types, and debris bins) is always assumed.

1) update_period
2) upper_N

### Safe to change (with percausions)
The following variables can be changed manually by the user, as long as certain conditions are met.

1) min_lifetime/CD/self.m0/min_dt/max_dt/dtfactor/t_max/setF107 : These can all be changed safely, but will not kick in until the next time NCell.update_lifetimes() is called. If you want to call it immediately, use the following code snippit (assuming the NCell instance is called system)
```Python3
system.update_lifetimes(system.t[system.time])
system.lupdate_time = system.time
```
2) sat_coll_probability_tables/rb_coll_probability_tables/sat_expl_probability_tables/rb_expl_probability_tables : These can techincally be updated at any time, as long as they contain valid probabilities, the cuumulative probability in each table is at most 1, and they have the shape (num_cells, num_cell, num_L, num_chi). However, it's heavily reccomended that you don't do this unless you familiar with the entirety of the code.

### Unsafe to change
Any variables not mentioned above could, in theory, be safely changed by the user, but this can only be garunteed if you're familiar with the entirety of the code.

### Accessing Cells
The cells contained in an NCell object are stored in the member list cells, ordered from lowest to highest altitude. They can be accessed directly using this parameter.

## Event Class
There's only one additional thing to cover in the Event class, the member variable last_event. This records the last time a frequency-based event occured (in yr), and is set to zero by default. The NCell class handles updating it while running simulations. If you want to introduce an Event into the class after running a simulation, you may need to set last_event to whatever time is right.

## But what if I want to...
There are many things that you could feasably want to change about the system which require changing the code itself. This is an overview of how to do a few of them (or why you really shouldn't)

### Change the Atmospheric Drag Lifetimes Model
Of all the possiblities outlined here, this is probably the easiest. The drag lifetime is always computed as the time it takes an object to decay from the top to the bottom of a cell. The function is called six times in the code (three times in NCell.\_\_init\_\_(), and three times in NCell.update_lifetimes()), and simply replacing these instances with your function should be sufficiant. If you want to change the currently implemented model, say by updating the solar cycle data, it's reccomended that you read the model documentation, and specifically look at the CIRA report it references.

### Turn Off Collisions
If you want to turn off a particular type of collision, the easiest way is to go into Cell.dxdt() and set the corresponding collision matrix to be all zeros. The collision matrices, and their shapes, are
1) dSdt/dDdt : Collisions between live/derelict satellites and debris, shape (num_sat, num_L, num_chi)
2) dSSdt/dSDdt/dDDdt : Collisions between satellites, where the first letter is the first satellite type and the second is the second satellite type, shape (num_sat, num_sat)
3) dSRdt/dDRdt : Collisions between live/derelict satellites and rocket bodies, shape (num_sat, num_rb)
4) dRdt : Collisions between rocket bodies and debris, shape (num_rb, num_L, num_chi)
5) dRRdt : Collisions between rocket bodies, shape (num_rb, num_rb)

### Turn Off Satellites/Rockets De-Orbiting and Decaying
This cannot be done in the current versions of the code. See the other implementations of simpler models.

### Change the Way Probabilities are Calculated
The entire model was designed around the NASA standard breakup model. Don't to this.

### Account for very High/Low Characteristic Lengths and A/M Ratios
This can be done, but you'll need to up the number of directions used in the probability table calculation by a lot. For A/M ratios outside the given range, it can take num_dir >= 100000 for accurate results to be computed. Essentially, it'll take a while. Also, due to a shortcut the calculation takes, there's a non-zero chance that your computer may run out of memory if num_dir is high enough.