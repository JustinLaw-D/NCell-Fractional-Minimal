# implementation of the NASA standard breakup model

from scipy.special import erf, erfinv
from scipy.stats import poisson
import numpy as np

def rand_poisson(ave, mx=np.inf):
    '''
    generates a random number from a poisson distribution (in the magnitude of the number), 
    retrying if the number is larger than mx

    Parameter(s):
    ave : expectation value of the Poisson distribution (can be negative)

    Keyword Parameter(s):
    mx : maximum desired random value (must be positive)

    Output(s):
    num : random value from poisson distribution
    '''

    if ave == 0 : return ave # nothing to do in this case
    sign_fac = round(ave/ave) # factor to account for the sign of the number
    first = True
    num = 0
    if abs(ave) > mx: 
        if mx > 10: # for small values things are probably fine
            print('WARNING: Time step may be too large')
        return mx*sign_fac
    while num > mx or first: # make sure the number isn't too large
        if first : first = False
        num = poisson.rvs(abs(ave))
    return num*sign_fac

def is_catastrophic(m_t, L, AM, v):
    '''
    determines if a collision between debris and a target is catastrophic

    Parameter(s):
    m_t : mass of the target (kg)
    L : characteristic length of the debris (m)
    AM : area to mass ratio of the debris (m^2/kg)
    v : relative velocity of the objects (km/s)

    Keyword Parameter(s): None

    Output(s):
    cat : True if the collision is catastrophic, False otherwise
    '''
    
    v *= 1000 # convert to m/s
    m_d = find_A(L)/AM # mass of the debris (on average)
    k_d = 0.5*m_d*(v**2) # relative kinetic energy of the debris
    dec_fact = (k_d/m_t)/1000 # factor for making the decision (J/g)
    if dec_fact >= 40 : return True
    else : return False

def calc_M(m_t, m_d, v):
    '''
    calculates the M factor used for L distribution calculation
    
    Parameter(s):
    m_t : target mass (array, kg)
    m_d : mass of the debris (array or 2d array, kg)
    v : collision velocity (km/s)

    Keyword parameter(s): None

    Output(s):
    M : value of M parameter for each combination (variable units)

    Note(s): output dimension order is (m_t dim, m_d dim)
    '''

    mt_size = m_t.size
    if (len(m_d.shape) == 1):
        m_t = np.resize(m_t, (m_d.size, mt_size)).transpose()
        m_d = np.resize(m_d, (mt_size, m_d.size))
    elif (len(m_d.shape) == 2):
        m_t = np.moveaxis(np.resize(m_t, (m_d.shape[1], m_d.shape[0], mt_size)), 2, 0)
        m_d = np.resize(m_d, (mt_size, m_d.shape[0], m_d.shape[1]))

    E_p = (0.5*m_d*((v*1000)**2)/m_t)/1000 # E_p in J/g
    M = np.empty(E_p.shape)
    M[E_p >= 40] = m_t[E_p >= 40] + m_d[E_p >= 40] # catestrophic collision
    M[E_p < 40] = m_d[E_p < 40]*v # non-catestrophic collision
    return M

def calc_Ntot(M, Lmin, Lmax, typ, C=1):
    '''
    calculates the total number of debris produced with characteristic length
    between Lmin and Lmax

    Parameter(s):
    M : fit parameter given by calc_M, ignored for explosions (variable units)
    Lmin : minimum characteristic length (m)
    Lmax : maximum characteristic length (m)
    typ : one of 'coll' (collision) or 'expl' (explosion)

    Keyword Parameter(s):
    C : fit parameter for explosions, ignored for collisions (default 1)

    Output(s):
    N : total number of fragments of Lmax > size > Lmin

    Note(s): returns 0 on an invalid type, M, C can be an arrays of any dimension
    '''

    if typ == 'coll' : return 0.1*(M**0.75)*((Lmin**(-1.71)) - (Lmax**(-1.71)))
    elif typ == 'expl' : return 6*C*(Lmin**(-1.6) - Lmax**(-1.6))
    else:
        print('WARNING: Invalid Debris Generation Type')
        return 0

def find_A(L):
    '''
    calculates the average cross-sectional area of an object with a given characteristic length

    Parameter(s):
    L : characteristic length (m)

    Keyword Parameter(s): None

    Output(s):
    A : average cross-sectional area (m^2)
    '''

    if L < 0.00167 : return 0.540424*(L**2)
    else : return 0.556945*(L**2.0047077)

def L_cdf(L, L_min, L_max, typ):
    '''
    calculates the cumulative distribution function for characteristic lengths
    from a collision/explosion at length L, assuming the distribution is truncated 
    at L_min and L_max

    Parameter(s):
    L : characteristic length (m)
    L_min : minimum characteristic length to consider (m)
    L_max : maximum characteristic length to consider (m)
    typ : one of 'coll' (collision) or 'expl' (explosion)

    Keyword Parameter(s): None

    Output(s):
    P : value of CDF at L

    Note(s): returns 0 on an invalid type, can take array for L
    '''

    if typ == 'coll' : beta = -1.71
    elif typ == 'expl' : beta = -1.6
    else:
        print('WARNING: Invalid Debris Generation Type')
        return 0
    return (L_min**beta - L**beta)/(L_min**beta - L_max**beta)

def X_cdf(x, x_min, x_max, L, typ):
    '''
    calculates the cumulative distribution function for log10(A/M) from a
    collision/explosion at values x and lengths L, assuming the distribution
    is truncated at x_min and x_max

    Parameter(s):
    x : log10(A/M) (array, log10(m^2/kg))
    x_min : minimum log10(A/M) value to consider (log10(m^2/kg))
    x_max : maximum log10(A/M) value to consider (log10(m^2/kg))
    L : characteristic length of the debris (array, m)
    typ : one of 'sat' (satellite) or 'rb' (rocket body)

    Keyword Parameter(s): None

    Output(s):
    P : value of CDF at x, L (first dim L, second x)

    Note(s): returns 0 on an invalid type, must take arrays for x, L
    '''

    if typ != 'sat' and typ != 'rb':
        print('WARNING: Invalid Debris Generator Type')
        return 0
    filter_top = (L >= 11/100) # get filters for which function to use
    filter_bottom = (L <= 8/100)
    filter_middle = ((filter_top | filter_bottom) == False)
    result = np.zeros((len(L), len(x)), dtype=np.double) # final result
    result[filter_top, :] = _X_cdf_11(x, x_min, x_max, L[filter_top], typ)
    result[filter_bottom, :] = _X_cdf_8(x, x_min, x_max, L[filter_bottom])
    lam_min, lam_max = np.log10(8/100), np.log10(11/100)
    P = (np.log10(L[filter_middle])-lam_min)/(lam_max-lam_min)
    result[filter_middle, :] = P*_X_cdf_11(x, x_min, x_max, L[filter_middle], typ) + (1-P)*_X_cdf_8(x, x_min, x_max, L[filter_middle])
    return result

def _X_cdf_8(x, x_min, x_max, L):
    '''
    calculates the cumulative distribution function for log10(A/M) from both collisions
    and explosions at values x and lengths L<=8cm, assuming the distribution is truncated 
    at x_min and x_max

    Parameter(s):
    x : log10(A/M) (array, log10(m^2/kg))
    x_min : minimum log10(A/M) value to consider (log10(m^2/kg))
    x_max : maximum log10(A/M) value to consider (log10(m^2/kg))
    L : characteristic length of the debris (array, m)

    Keyword Parameter(s): None

    Output(s):
    P : values of CDF at x, L (first dim L, second x)
    '''
    lam = np.log10(L)

    # define functions for determining normal distribution parameters
    def mu_soc(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        bot_filer = (lambda_c <= -1.75)
        mid_filter = (lambda_c > -1.75) & (lambda_c < -1.25)
        top_filter = lambda_c >= -1.25
        result[bot_filer] = -0.3
        result[mid_filter] = -0.3 - 1.4*(lambda_c[mid_filter] + 1.75)
        result[top_filter] = -1
        return result

    def sigma_soc(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        top_filter = (lambda_c > -3.5)
        result[top_filter == False] = 0.2
        result[top_filter] = 0.2 + 0.1333*(lambda_c[top_filter] + 3.5)
        return result
    
    mu = mu_soc(lam) # calculate parameters
    sigma = sigma_soc(lam)
    x_2d = np.resize(x, (len(L), len(x)))
    x_2d = x_2d.transpose() # needs to be in this direction to broadcast correctly in the next calculations
    C = 1/(erf((x_max-mu)/(np.sqrt(2)*sigma)) - erf((x_min-mu)/(np.sqrt(2)*sigma))) # normalization factors for each L
    # compute total distribution
    return np.swapaxes(C*(erf((x_2d-mu)/(np.sqrt(2)*sigma)) - erf((x_min-mu)/(np.sqrt(2)*sigma))), 0, 1)

def _X_cdf_11(x, x_min, x_max, L, typ):
    '''
    calculates the cumulative distribution function for log10(A/M) from collisions/explosions 
    at values x and lengths L>=11cm, assuming the distribution is truncated at x_min 
    and x_max

    Parameter(s):
    x : log10(A/M) (array, log10(m^2/kg))
    x_min : minimum log10(A/M) value to consider (log10(m^2/kg))
    x_max : maximum log10(A/M) value to consider (log10(m^2/kg))
    L : characteristic length of the debris (array, m)
    typ : one of 'sat' (satellite) or 'rb' (rocket body)

    Keyword Parameter(s): None

    Output(s):
    P : value of CDF at x, L (first dim L, second x)

    Note(s): returns 0 on an invalid type
    '''
    lam = np.log10(L)

    # define functions for determining normal distribution parameters
    def alpha_sc(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        bot_filter = (lambda_c <= -1.95)
        mid_filter = (lambda_c > -1.95) & (lambda_c < 0.55)
        top_filter = (lambda_c >= 0.55)
        result[bot_filter] = 0
        result[mid_filter] = 0.3 + 0.4*(lambda_c[mid_filter] + 1.2)
        result[top_filter] = 1
        return result

    def mu1_sc(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        bot_filter = (lambda_c <= -1.1)
        mid_filter = (lambda_c > -1.1) & (lambda_c < 0)
        top_filter = (lambda_c >= 0)
        result[bot_filter] = -0.6
        result[mid_filter] = -0.6 - 0.318*(lambda_c[mid_filter] + 1.1)
        result[top_filter] = -0.95
        return result

    def sigma1_sc(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        bot_filter = (lambda_c <= -1.3)
        mid_filter = (lambda_c > -1.3) & (lambda_c < -0.3)
        top_filter = (lambda_c >= -0.3)
        result[bot_filter] = 0.1
        result[mid_filter] = 0.1 + 0.2*(lambda_c[mid_filter] + 1.3)
        result[top_filter] = 0.3
        return result

    def mu2_sc(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        bot_filter = (lambda_c <= -0.7)
        mid_filter = (lambda_c > -0.7) & (lambda_c < -0.1)
        top_filter = (lambda_c >= -0.1)
        result[bot_filter] = -1.2
        result[mid_filter] = -1.2 - 1.333*(lambda_c[mid_filter] + 0.7)
        result[top_filter] = -2
        return result

    def sigma2_sc(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        bot_filter = (lambda_c <= -0.5)
        mid_filter = (lambda_c > -0.5) & (lambda_c < -0.3)
        top_filter = (lambda_c >= -0.3)
        result[bot_filter] = 0.5
        result[mid_filter] = 0.5 - (lambda_c[mid_filter] + 0.5)
        result[top_filter] = 0.3
        return result

    def alpha_rb(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        bot_filter = (lambda_c <= -1.4)
        mid_filter = (lambda_c > -1.4) & (lambda_c < 0)
        top_filter = (lambda_c >= 0)
        result[bot_filter] = 1
        result[mid_filter] = 1 - 0.3571*(lambda_c[mid_filter] + 1.4)
        result[top_filter] = 0.5
        return result

    def mu1_rb(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        bot_filter = (lambda_c <= -0.5)
        mid_filter = (lambda_c > -0.5) & (lambda_c < 0)
        top_filter = (lambda_c >= 0)
        result[bot_filter] = -0.45
        result[mid_filter] = -0.45 - 0.9*(lambda_c[mid_filter] + 0.5)
        result[top_filter] = -0.9
        return result

    def sigma1_rb(lambda_c):
        return np.full(len(lambda_c), 0.55)

    def mu2_rb(lambda_c):
        return np.full(len(lambda_c), -0.9)

    def sigma2_rb(lambda_c):
        result = np.empty(len(lambda_c), dtype=np.double)
        bot_filter = (lambda_c <= -1)
        mid_filter = (lambda_c > -1) & (lambda_c < 0.1)
        top_filter = (lambda_c >= 0.1)
        result[bot_filter] = 0.28
        result[mid_filter] = 0.28 - 0.1636*(lambda_c[mid_filter] + 1)
        result[top_filter] = 0.1
        return result
    
    if typ == 'sat':
        mu1 = mu1_sc(lam) # calculate parameters
        sigma1 = sigma1_sc(lam)
        mu2 = mu2_sc(lam)
        sigma2 = sigma2_sc(lam)
        alpha = alpha_sc(lam)
    elif typ == 'rb':
        mu1 = mu1_rb(lam) # calculate parameters
        sigma1 = sigma1_rb(lam)
        mu2 = mu2_rb(lam)
        sigma2 = sigma2_rb(lam)
        alpha = alpha_rb(lam)
    else:
        print('WARNING: Invalid Debris Generator Type')
        return 0
    x_2d = np.resize(x, (len(L), len(x)))
    x_2d = x_2d.transpose() # needs to be in this direction to broadcast correctly in the next calculations
    # compute normalization factors for each L
    top = alpha*erf((x_max-mu1)/(np.sqrt(2)*sigma1)) + (1-alpha)*erf((x_max-mu2)/(np.sqrt(2)*sigma2))
    bot = alpha*erf((x_min-mu1)/(np.sqrt(2)*sigma1)) + (1-alpha)*erf((x_min-mu2)/(np.sqrt(2)*sigma2))
    C = 1/(top - bot)
    # compute total distribution
    fac_one = erf((x_2d-mu1)/(np.sqrt(2)*sigma1)) - erf((x_min-mu1)/(np.sqrt(2)*sigma1))
    fac_two = erf((x_2d-mu2)/(np.sqrt(2)*sigma2)) - erf((x_min-mu2)/(np.sqrt(2)*sigma2))
    return np.swapaxes(C*(alpha*fac_one + (1-alpha)*fac_two),0,1)

def v_cdf(v, x, typ):
    '''
    evaluates cdf for log10(Delta v) values at given vs and xs

    Parameter(s):
    v : log10(delta v) value to evaluate at (array, log10(m/s))
    x : log10(A/M) value of the debris (array, log10(m^2/kg))
    typ : one of 'coll' (collision) or 'expl' (explosion)

    Keyword Parameter(s): None

    Output(s):
    P : value of the CDF at xs, vs (first dim x, second v)

    Note(s): returns 0 on an invalid type
    '''
    
    if typ == 'coll' : mu = 0.9*x + 2.9 # calculate normal distribution parameters
    elif typ == 'expl' : mu = 0.2*x + 1.85
    else:
        print('WARNING: Invalid Debris Generation Type')
        return 0
    sigma_fac = 0.4*np.sqrt(2)
    C = 1/2 # calculate normalization factor
    v_2d = np.resize(v, (len(x), len(v))) # 2d array of v
    v_2d = v_2d.transpose() # needed for proper broadcasting
    # calculate CDF value
    return np.swapaxes(C*(erf((v_2d-mu)/sigma_fac) + 1),0,1)

def vprime_cdf(V, v0, theta, phi, x, typ):
    '''
    evaluates cdf for the post-collision speed V, given a pre-collision
    orbital speed v0, post-collision directions of (theta, phi), and xs

    Parameter(s):
    V : post-collison speed to evaluate at (m/s)
    v0 : pre-collision orbital speed (m/s)
    theta : inclination angle (array, rad)
    phi : azimuthal angle (array, rad)
    x : log10(A/M) value of the debris (array, log10(m^2/kg))
    typ : one of 'coll' (collision) or 'expl' (explosion)

    Keyword Parameter(s): None

    Output(s):
    P : value of the CDF at V (first dim x, second direction)

    Note(s): returns 0 on an invalid type
    '''
    result = np.zeros((len(x), len(theta)), dtype=np.double)
    descriminate = (v0*np.sin(theta)*np.cos(phi))**2 - (v0**2-V**2)
    # descriminate < 0 means a probability of zero
    del_v_max = -v0*np.sin(theta)*np.cos(phi) + np.sqrt(descriminate)
    del_v_min = -v0*np.sin(theta)*np.cos(phi) - np.sqrt(descriminate)
    # del_v_max < 0 is also not possible
    possible = (np.isnan(del_v_max) == False) & (del_v_max >=0)
    del_v_min_low = possible & (del_v_min <= 0)
    del_v_min_high = possible & (del_v_min_low == False)
    result[:, del_v_min_low] = v_cdf(np.log10(del_v_max[del_v_min_low]), x, typ)
    result[:, del_v_min_high] = v_cdf(np.log10(del_v_max[del_v_min_high]), x, typ) - v_cdf(np.log10(del_v_min[del_v_min_high]), x, typ)
    return result

# ALL FUNCTIONS PAST THIS POINT ARE NO LONGER USED, AND MAY OR MAY NOT BE FUNCTIONAL

def randL(num, L_min, L_max, typ):
    '''
    generates num random characteristic lengths for debris from a collision/explosion

    Parameter(s):
    num : number of random lengths to generate
    L_min : minimum characteristic length to consider (m)
    L_max : maximum characteristic length to consider (m)
    typ : one of 'coll' (collision) or 'expl' (explosion)

    Keyword Parameter(s): None

    Output(s):
    L : array of random characteristic lengths (m)

    Note(s): returns 0 on an invalid type
    '''

    lam_min, lam_max = np.log10(L_min), np.log10(L_max)
    if typ == 'coll' : beta = -1.71
    elif typ == 'expl' : beta = -1.6
    else:
        print('WARNING: Invalid Debris Generation Type')
        return 0
    P = np.random.uniform(size=num) # get random P values
    lam = np.log10(10**(beta*lam_min) - P*(10**(beta*lam_min) - 10**(beta*lam_max)))/beta
    return 10**lam

def randX(num, x_min, x_max, L, typ):
    '''
    generates num random log10(A/M) values for debris from a collision/explosion

    Parameter(s):
    num : number of random values to generate
    x_min : minimum log10(A/M) value to consider (log10(m^2/kg))
    x_max : maximum log10(A/M) value to consider (log10(m^2/kg))
    L : characteristic length of the debris (m)
    typ : one of 'sat' (satellite) or 'rb' (rocket body)

    Keyword Parameter(s): None

    Output(s):
    x : array of random log10(A/M) values (log10(m^2/kg))

    Note(s): returns 0 on an invalid type
    '''

    if typ != 'sat' and typ != 'rb':
        print('WARNING: Invalid Debris Generator Type')
        return 0
    if L >= 11/100 : return _randX_11(num, x_min, x_max, L, typ)
    elif L <= 8/100 : return _randX_8(num, x_min, x_max, L)
    else:
        if typ == 'sat' : comp = 10*(np.log10(L) + 1.05)
        else : comp = 10*(np.log10(L) + 1.76)
        if np.random.uniform() > comp : return _randX_11(num, x_min, x_max, L, typ)
        else : return _randX_8(num, x_min, x_max, L)

def _randX_8(num, x_min, x_max, L):
    '''
    generates num random log10(A/M) values for debris from a collision/explosion, 
    assuming that the characteristic length of the debris is less than 8cm

    Parameter(s):
    num : number of random values to generate
    x_min : minimum log10(A/M) value to consider (log10(m^2/kg))
    x_max : maximum log10(A/M) value to consider (log10(m^2/kg))
    L : characteristic length of the debris (m)

    Keyword Parameter(s): None

    Output(s):
    x : array of random log10(A/M) values (log10(m^2/kg))
    '''

    lam = np.log10(L)

    # define functions for determining normal distribution parameters
    def mu_soc(lambda_c):
        if lambda_c <= -1.75 : return -0.3
        elif lambda_c < -1.25 : return -0.3 - 1.4*(lambda_c + 1.75)
        else : return -1

    def sigma_soc(lambda_c):
        if lambda_c <= -3.5 : return 0.2
        else : return 0.2 + 0.1333*(lambda_c + 3.5)
    
    mu = mu_soc(lam) # calculate parameters
    sigma = sigma_soc(lam)
    C = 1/(erf((x_max-mu)/(np.sqrt(2)*sigma)) - erf((x_min-mu)/(np.sqrt(2)*sigma))) # normalization factor
    P = np.random.uniform(size=num) # get random P values
    # use these to generate random x-values
    x = sigma*np.sqrt(2)*erfinv(P/C + erf((x_min - mu)/(sigma*np.sqrt(2)))) + mu
    return x

def _randX_11(num, x_min, x_max, L, typ):
    '''
    generates num random log10(A/M) values for debris from a collision/explosion, 
    assuming that the characteristic length of the debris is greater than 11cm

    Parameter(s):
    num : number of random values to generate
    x_min : minimum log10(A/M) value to consider (log10(m^2/kg))
    x_max : maximum log10(A/M) value to consider (log10(m^2/kg))
    L : characteristic length of the debris (m)
    typ : one of 'sat' (satellite) or 'rb' (rocket body)

    Keyword Parameter(s): None

    Output(s):
    x : array of random log10(A/M) values (log10(m^2/kg))

    Note(s): returns 0 on an invalid type
    '''

    lam = np.log10(L)

    # define functions for determining normal distribution parameters
    def alpha_sc(lambda_c):
            if lambda_c <= -1.95 : return 0
            elif lambda_c < 0.55 : return 0.3 + 0.4*(lambda_c + 1.2)
            else : return 1

    def mu1_sc(lambda_c):
        if lambda_c <= -1.1 : return -0.6
        elif lambda_c < 0 : return -0.6 - 0.318*(lambda_c + 1.1)
        else : return -0.95

    def sigma1_sc(lambda_c):
        if lambda_c <= -1.3 : return 0.1
        elif lambda_c < -0.3 : return 0.1 + 0.2*(lambda_c + 1.3)
        else : return 0.3

    def mu2_sc(lambda_c):
        if lambda_c <= -0.7 : return -1.2
        elif lambda_c < -0.1 : return -1.2 - 1.333*(lambda_c + 0.7)
        else : return -2

    def sigma2_sc(lambda_c):
        if lambda_c <= -0.5 : return 0.5
        elif lambda_c < -0.3 : return 0.5 - (lambda_c + 0.5)
        else : return 0.3

    def alpha_rb(lambda_c):
            if lambda_c <= -1.4 : return 1
            elif lambda_c < 0 : return 1 - 0.3571*(lambda_c + 1.4)
            else : return 0.5

    def mu1_rb(lambda_c):
        if lambda_c <= -0.5 : return -0.45
        elif lambda_c < 0 : return -0.45 - 0.9*(lambda_c + 0.5)
        else : return -0.9

    def sigma1_rb(lambda_c):
        return 0.55

    def mu2_rb(lambda_c):
        return -0.9

    def sigma2_rb(lambda_c):
        if lambda_c <= -1 : return 0.28
        elif lambda_c < 0.1 : return 0.28 - 0.1636*(lambda_c + 1)
        else : return 0.1
    
    if typ == 'sat':
        mu1 = mu1_sc(lam) # calculate parameters
        sigma1 = sigma1_sc(lam)
        mu2 = mu2_sc(lam)
        sigma2 = sigma2_sc(lam)
        alpha = alpha_sc(lam)
    elif typ == 'rb':
        mu1 = mu1_rb(lam) # calculate parameters
        sigma1 = sigma1_rb(lam)
        mu2 = mu2_rb(lam)
        sigma2 = sigma2_rb(lam)
        alpha = alpha_rb(lam)
    else:
        print('WARNING: Invalid Debris Generator Type')
        return 0
    # compute normalization factor
    top = alpha*erf((x_max-mu1)/(np.sqrt(2)*sigma1)) + (1-alpha)*erf((x_max-mu2)/(np.sqrt(2)*sigma2))
    bot = alpha*erf((x_min-mu1)/(np.sqrt(2)*sigma1)) + (1-alpha)*erf((x_min-mu2)/(np.sqrt(2)*sigma2))
    C = 1/(top - bot)
    x_table = np.linspace(x_min, x_max, num=1000) # table of x values
    # corresponding table of P values
    P_table = C*(alpha*erf((x_table-mu1)/(np.sqrt(2)*sigma1)) + (1-alpha)*erf((x_table-mu2)/(np.sqrt(2)*sigma2)) - bot)
    P = np.random.uniform(size=num) # get random P values
    # use these to generate random x-values
    x = np.zeros(P.shape)
    for i in range(len(P)):
        index = np.abs(P_table - P[i]).argmin() # find location of closest value on the table
        x[i] = x_table[index] # use this to find the corresponding x-value
    return x

def randv(num, x, typ):
    '''
    generates num random log10(Delta v) values for debris from a collision/explosion

    Parameter(s):
    num : number of random values to generate
    x : log10(A/M) value of the debris (log10(m^2/kg))
    typ : one of 'coll' (collision) or 'expl' (explosion)

    Keyword Parameter(s): None

    Output(s):
    v : array of random log10(Delta v) values (log10(m/s))

    Note(s): returns 0 on an invalid type
    '''

    if typ == 'coll' : mu = 0.9*x + 2.9 # calculate normal distribution parameters
    elif typ == 'expl' : mu = 0.2*x + 1.85
    else:
        print('WARNING: Invalid Debris Generation Type')
        return 0
    sigma_fac = 0.4*np.sqrt(2)
    C = 1/2 # calculate normalization factor
    P = np.random.uniform(size=num) # get random P values
    # use these to generate random v-values
    v = sigma_fac*erfinv(P/C - 1) + mu
    return v

def rand_direction(num):
    '''
    generates random directions in 3-d space

    Parameter(s):
    num : number of random values to generate

    Keyword Parameter(s): None

    Output(s)
    u : array of 3-d unit vectors in cartesian coordinates

    Notes: Output u is a matrix, where the first row a list of x-coordinates,
    the second y-coordinates, and the third z-coordinates
    '''

    theta = np.random.uniform(0.0, np.pi, size=num) # generate inclinations
    phi = np.random.uniform(0.0, 2*np.pi, size=num) # generate azimuthal angles
    to_return = np.zeros((3, num)) # first row is x, second is y, third is z
    to_return[0, :] = np.cos(phi)*np.sin(theta)
    to_return[1, :] = np.sin(phi)*np.sin(theta)
    to_return[2, :] = np.cos(theta)
    return to_return