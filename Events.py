# classes for discrete events

import numpy as np

class Event:
    
    def __init__(self, alt, time=None, freq=None):
        '''
        constructor for general event class

        Paremeter(s):
        alt : altitude of the event (km)

        Keyword Parameter(s):
        time : list of times that the event occurs (yr, default None)
        freq : frequency the event occurs at (1/yr, default None)

        Output(s): instance of Event

        Note(s): time and freq cannot both be None
        '''

        if time is None and freq is None:
            print('Invlid Event : No occurance time specified')

        self.time = time
        self.last_event = 0 # time of the last event (yr)
        self.freq = freq
        self.alt = alt

    def run_event(self, S, S_d, D, R, N, logL_edges, chi_edges):
        '''
        function representing the discrete event occuring

        Input(s):
        S : number of live satellites of each type in the current cell (list of floats)
        S_d : number of de-orbiting satellites of each type in the current cell (list of floats)
        D : number of derelict satellites of each type in the current cell (list of floats)
        R : number of rocket bodies of each type in the current cell (list of floats)
        N : binned amount of debris in current cell (2d array)
        logL_edges : logL edge values for the bins (log10(m))
        chi_edges : chi edge values for the bins (log10(m^2/kg))

        Keyword Input(s): None

        Output(s):
        dS : change in the number of live satellites of each type in the current cell (array of floats)
        dS_d : change in the number of de-orbiting satellites of each type in the current cell (array of floats)
        dD : change in the number of derelict satellites of each type in the cell (array of floats)
        dR : change in the number of rocket bodies of each type in the cell (array of floats)
        dN : change in the number of debris in the curren cell, not including debris
             produced by collisions
        coll : list of collisions occuring in the current cell in the form [(kg, kg, typ, #)],
               i.e. [(m1, m2, typ, number of collisions)]. typ can be one of 'sat' (satellite-satellite),
               'sr' (satellite-rocket, where satellite is m1), or 'rb' (rocket-rocket)
        expl : list of explosions occuring in the current cell in the form [(C, typ, #)], where
               C is the relevant fit constant and typ is the type of body exploding ('sat' or 'rb)

        Note(s): this function is meant to be overwritten, and in the default form just returns
                 zero
        '''

        return np.zeros(S.shape), np.zeros(S_d.shape), np.zeros(D.shape), np.zeros(R.shape), np.zeros(N.shape), [], []

# class for handling basic explosions
class ExplEvent(Event):
    
    def __init__(self, alt, expl_list, time=None, freq=None):
        '''
        constructor for a basic explosions event class

        Paremeter(s):
        alt : altitude of the event (km)
        expl_list : list of explosions occuring in the current cell on an event in the form [(C, typ, #)], 
                    where C is the relevant fit constant and typ is the type of body exploding ('sat' or 'rb)

        Keyword Parameter(s):
        time : list of times that the event occurs (yr, default None)
        freq : frequency the event occurs at (1/yr, default None)

        Output(s): instance of Event

        Note(s): time and freq cannot both be None
        '''

        super().__init__(alt, time=time, freq=freq)
        self.expl_list = expl_list

    def run_event(self, S, S_d, D, R, N, logL_edges, chi_edges):
        '''
        function representing the discrete event occuring

        Input(s):
        S : number of live satellites of each type in the current cell (list of floats)
        S_d : number of de-orbiting satellites of each type in the current cell (list of floats)
        D : number of derelict satellites of each type in the current cell (list of floats)
        R : number of rocket bodies of each type in the current cell (list of floats)
        N : binned amount of debris in current cell (2d array)
        logL_edges : logL edge values for the bins (log10(m))
        chi_edges : chi edge values for the bins (log10(m^2/kg))

        Keyword Input(s): None

        Output(s):
        dS : change in the number of live satellites of each type in the current cell (array of floats)
        dS_d : change in the number of de-orbiting satellites of each type in the current cell (array of floats)
        dD : change in the number of derelict satellites of each type in the cell (array of floats)
        dR : change in the number of rocket bodies of each type in the cell (array of floats)
        dN : change in the number of debris in the curren cell, not including debris
             produced by collisions
        coll : list of collisions occuring in the current cell in the form [(kg, kg, typ, #)],
               i.e. [(m1, m2, typ, number of collisions)]. typ can be one of 'sat' (satellite-satellite),
               'sr' (satellite-rocket, where satellite is m1), or 'rb' (rocket-rocket)
        expl : list of explosions occuring in the current cell in the form [(C, typ, #)], where
               C is the relevant fit constant and typ is the type of body exploding ('sat' or 'rb)

        Note(s): this function is meant to be overwritten, and in the default form just returns
                 zero
        '''

        return np.zeros(S.shape), np.zeros(S_d.shape), np.zeros(D.shape), np.zeros(R.shape), np.zeros(N.shape), [], self.expl_list

# class for handling basic collisions
class CollEvent(Event):
    
    def __init__(self, alt, coll_list, time=None, freq=None):
        '''
        constructor for a basic collisions event class

        Paremeter(s):
        alt : altitude of the event (km)
        coll_list : list of collisions occuring in the current cell in the form [(kg, kg, typ, #)],
                    i.e. [(m1, m2, typ, number of collisions)]. typ can be one of 'sat' (satellite-satellite),
                    'sr' (satellite-rocket, where satellite is m1), or 'rb' (rocket-rocket)

        Keyword Parameter(s):
        time : list of times that the event occurs (yr, default None)
        freq : frequency the event occurs at (1/yr, default None)

        Output(s): instance of Event

        Note(s): time and freq cannot both be None
        '''

        super().__init__(alt, time=time, freq=freq)
        self.coll_list = coll_list

    def run_event(self, S, S_d, D, R, N, logL_edges, chi_edges):
        '''
        function representing the discrete event occuring

        Input(s):
        S : number of live satellites of each type in the current cell (list of floats)
        S_d : number of de-orbiting satellites of each type in the current cell (list of floats)
        D : number of derelict satellites of each type in the current cell (list of floats)
        R : number of rocket bodies of each type in the current cell (list of floats)
        N : binned amount of debris in current cell (2d array)
        logL_edges : logL edge values for the bins (log10(m))
        chi_edges : chi edge values for the bins (log10(m^2/kg))

        Keyword Input(s): None

        Output(s):
        dS : change in the number of live satellites in the current cell (array of floats)
        dS_d : change in the number of de-orbiting satellites in the current cell (array of floats)
        dD : change in the number of derelict satellites in the cell (array of floats)
        dR : change in the number of rocket bodies of each type in the cell (array of floats)
        dN : change in the number of debris in the curren cell, not including debris
             produced by collisions
        coll : list of collisions occuring in the current cell in the form [(kg, kg, typ, #)],
               i.e. [(m1, m2, typ, number of collisions)]. typ can be one of 'sat' (satellite-satellite),
               'sr' (satellite-rocket, where satellite is m1), or 'rb' (rocket-rocket)
        expl : list of explosions occuring in the current cell in the form [(C, typ, #)], where
               C is the relevant fit constant and typ is the type of body exploding ('sat' or 'rb)

        Note(s): this function is meant to be overwritten, and in the default form just returns
                 zero
        '''

        return np.zeros(S.shape), np.zeros(S_d.shape), np.zeros(D.shape), np.zeros(R.shape), np.zeros(N.shape), self.coll_list, []
