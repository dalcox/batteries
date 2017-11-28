#analytically combined equation
def simplebattfunc(x, i):
    """
    Analytically solved current and voltage within electrode
    """
    i0, i1 = i
    di = i1
    d2i = a*io*(n*F)/(R*T)*(-I/s + i0*(1/s + 1/K))
    return di, d2i

#systems of differential equations with drop-in replacements possible
def linearbattfunc(x, IV):
    """
    Linear current and voltage within single electrode
    """
    i1, i2, V1, V2 = IV
    di2 = a*io*(n*F)/(R*T)*(V1 - V2)
    #Kinetics
    di1 = -di2
    #charge neutrality
    dV1 = -i1/s
    #solids ohms law
    dV2 = -i2/K
    #liquids ohms law
    return di1, di2, dV1, dV2

def BVbattfunc(x, IV):
    """
    Full Butler-Volmer kinetics
    """
    i1, i2, V1, V2 = IV
    di2 = a*io*(np.exp((n*F)/(R*T)*aa*(V1 - V2)) - np.exp((n*F)/(R*T)*-ac*(V1 - V2)))
    #Kinetics
    di1 = -di2
    #charge neutrality
    dV1 = -i1/s
    #solids ohms law
    dV2 = -i2/K
    #liquids ohms law
    return di1, di2, dV1, dV2

def Tafelfunc_a(x, IV):
    """
    uses Tafel kinetics
    """
    i1, i2, V1, V2 = IV

    taff = aa*n*F/(R*T)*(V1-V2)

    di2 = aa*io*np.exp(taff)
    #Kinetics
    di1 = -di2
    #charge neutrality
    dV1 = -i1/s
    #solids ohms law
    dV2 = -i2/K
    #liquids ohms law

    return di1, di2, dV1, dV2

def Tafelfunc_c(x, IV):
    """
    uses Tafel kinetics
    """
    i1, i2, V1, V2 = IV

    taff = -ac*n*F/(R*T)*(V1-V2)

    di2 = -a*io*np.exp(taff)
    #Kinetics
    di1 = -di2
    #charge neutrality
    dV1 = -i1/s
    #solids ohms law
    dV2 = -i2/K
    #liquids ohms law

    return di1, di2, dV1, dV2
