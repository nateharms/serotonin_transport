def multistep(c, vmax1, Km1, K1, vmax2, Km2, K2):
    """
    c is a vector of concentrations of Tryptophan, 5-hydroxytryptophoan, and Serotonin
    based on Michaelis Menten kinetics
    Returns the rate dependent on the concentrations
    """
    cTry, c5HTP, cSer = c
    k1fwd = vmax1/(Km1+cTry) #michaelus menton kinetics
    k1rev = k1fwd/K1 #equilibrium constant
    dcTry = -k1fwd*cTry + k1rev*c5HTP

    #Same reaction type for second reaction step
    k2fwd = vmax2/(Km2+c5HTP)
    k2rev = k2fwd/K2
    dcSer = k2fwd*c5HTP - k2rev*cSer

    dc5HTP = -dcTry - dcSer
    return dcTry, dc5HTP, dcSer
