import numpy as np
from reactionKinetics import multistep


class LaminarFlow:
    '''
    Model for laminar flow through intestines, including reactions
    '''
    def __init__(self, length, radius, max_velocity, serConditions, trypConditions, kinetics, rings, sections, time #M_Beta_file = 'valuesBetaM.csv'):
        """
        Initializes the variables to run begin running the simulation: Maybe we make this general so we can use the same class
        to run Serotonin and Tryptophan. Haven't really decided yet. Different compounds would change our Graetz number as well.
        Beta and M values come from the mathematical solution to the differential equations.

        Input Variables
        ------------------------------------------------------------------------------------
        length: the length of the intestines (m)
        radius: the radius of the intestines (m)
        max_velocity: velocity of the media (m/s)
        serConditions: type: NamedTuple, with ('Concentration' (mM), 'Diffusivity', 'Wall Permeability')
        trypConcentration: type: NamedTuple, with ('Concentration' (mM), 'Diffusivity', 'Wall Permeability')
        kinetics: (vmax1, Km1, K1, vmax2, Km2, K2) as variables for Michaelis Menten, hrs, mM
        rings: the number of ring partitions of intestines looked at (number of rings)
        sections: the number of sub-intervals of intestines looked at (number of sub-sections)
        M_Beta_file: the file that contains the M_Beta values needed for calculations
        """

        # Dimentions of the intestines
        self.length = length
        self.max_velocity = max_velocity
        self.dt = (length/max_velocity/3600)/n # to get hours
        self.radius = radius
        self.kinetics = kinetics.kinetics
        self.wallKinetics = kinetics.wallKinetics # help Jason
        
        # Subsection information
        self.rings = rings
        self.sections = sections
        
        # The initial concentrations for both of the tryp an ser.
        self.serConcentration = np.zeros(n)
        self.serConcentration[0] = serConditions.Concentration

        self.trypConcentration = np.zeros(n)
        self.trypConcentration[0] = trypConditions.Concentration

        self.HTPConcentration = np.zeros(n)

        # The diffusivities, wall permabilities, and effective permabilities for ser and tryp.
        self.serDiffusivity = serConditions.Diffusivity #6.2424e-8 m^2/sec
        self.serWallPerm = serConditions.Permeability #7.576e-13 m^2/sec
        self.serEffPerm = self.serWallPerm * self.radius / self.serDiffusivity

        self.trypDiffusivity = trypConditions.Diffusivity #5.386e-8 m^2/sec
        self.trypWallPerm = trypConditions.Permeability #6.44e-4 m^2/sec
        self.trypEffPerm = self.trypWallPerm * self.radius / self.trypDiffusivity
        
        self.htpDiffusivity = htpConditions.Diffusivity #4.995e-8 m^2/sec
        # self.htpWallPerm = htpConditions.Permeability
        # self.htpEffPerm = self.htpWallPerm * self.radius / self.htpDiffusivity

        # The MBeta values
        self.MBeta = np.genfromtxt(M_Beta_file, delimiter = ',', skip_header = 1)
        
        # self.getGraetz(n)
        self.getConcentration()

    def getGraetz(self, n):
        """
        This method is designed to determine the graetz number for both ser and tryp across the intestines.
        This value changes over the length of the intestines, however, not much else should.

        Inputs:
        ------------------------------------------------------------------------------------
        n: number of sub-intervals of the intestines

        Outputs:
        ------------------------------------------------------------------------------------
        lengths: a list of length n that has the legnths at which the Graetz numbers were calculated at
        serGraetz: a list of Graetz numbers for ser (length n)
        trypGraetz: a list of Graetz numbers for tryp (length n)
        """

        self.lengths = lengths = np.linspace(0,self.length,n)
        self.serGraetz = self.serDiffusivity/(self.max_velocity*self.radius**2)*lengths
        self.trypGraetz = self.serGraetz*self.trypDiffusivity/self.serDiffusivity
                 
    def crossSecAreas(self, rings, velocity_max):
        "gives an array of radial ring areas for a tube based on given radius and number of subsections"
        totarea=np.zeros(rings)
        velocityProfile=np.zeros(rings)
        r=self.radius/rings
        for m in range(rings):
            totarea[m]=np.pi*(r*(m+1))**2 - np.pi*(r*(m))**2
            velocityProfile[m] = velocity_max*(self.radius**2-(r*m)**2)    
        return totarea, velocityProfile

    def getConcentration(self):
        # Creating three 3-D arrays to keeps track of the following:
        # ~ Substance Concentration
        # ~ Ring number (in the r direction)
        # ~ Section number (in the z direction)

        # Conventions that are used throughout the script, named for convienence.
        rings = self.rings
        sections = self.sections
        time = self.time

        # 3-D Arrays
        trypconc = np.zeros((len(time), rings, sections))
        serotonconc = np.zeros_like(trypconc)
        htpconc = np.zeros_like(trypconc)

        # Adding the initial conditions to the tryp and ser arrays
        trypconc[0,:,0] = trypConditions.concentration
        serconc[0,:,0] = serConditions.concentration

        # Initializing the total amount of ser taken up.
        totalSerotoninUptake = 0

        # For loop to look through each time.
        for i in range(len(time)-1):
            #--------------------------------------------------------------------------------#
            # Diffusion Calculations

            # Determining the laplacians and first deravitives for each concentration curve.
            # Laplicians are used for diffusion calculations
            # First deravitives are used for convection calculations

            delZArrays = []
            delRArrays = []
            lapZArrays = []
            lapRArrays = []

            concentrationTuple = (trypconc[i,:,:], serconc[i,:,:], htpconc[i,:,:])
            for j, Z in enumerate(concentrationtuple): #0 - tryp, 1 - ser, 2 - htp
                delZArrays.append(deltaZ(Z, dz))
                delRArrays.append(deltaR(Z, r))
                lapZArrays.append(laplacianZ(Z, dz))
                lapRArrays.append(laplacianR(Z, r))

            #--------------------------------------------------------------------------------#
            # Reaction Calculations 

            # Conversion of tryp via reaction.
            # This is split into two parts, the reaction in the bulk and the reaction at the intestine wall.

            # Initialize the reaction arrays
            try_rxn = np.zeros_like(trypconc[i,:,:])
            ser_rxn = np.zeros_like(trypconc[i,:,:])
            htp_rxn = np.zeros_like(trypconc[i,:,:])


            for q in range(rings):
                if q == 0: # Reaction on the surface
                    try_rxn[q,:], ser_rxn[q,:], htp_rxn[q,:] = reactionKinetics.multistep(concentrationTuple[:][i,q,:], *self.wallKinetics) #Need to make a wallKinetics portion
                else: # Reaction in the bulk
                    try_rxn[q,:], ser_rxn[q,:], htp_rxn[q,:] = reactionKinetics.multistep(concentrationTuple[:][i,q,:], *self.kinetics)
            """
            #tryp_rxn, ser_rxn, htp_rxn = reactionKinetics.multistep(concentrationTuple, *self.kinetics) #trypconc, serotonconc, htpconc

            # Rxn at the surface
            #for Z in (tryp_rxn, ser_rxn, htp_rxn):
            #    Z[-1,:] *= 0 # Removing the existing reaction rates because they aren't correct for the surface.

            #for Z in concentrationTuple[0]:
            #    tryp_rxn[-1,:], ser_rxn[-1,:], htp_rxn[-1,:] = reactionKinetics.multistep(concentrationTuple)
            """
            #--------------------------------------------------------------------------------#
            #Boundary Condition Layer

            #### I do not know what is going on here ####
            for j, Z in enumerate(trypconc, serotonconc, htpconc): 
                Z[i,0,:] = Z[i,1,:]
                Z[i,-1,:] = 0 #insert some magical boundary condition here (wall permeability * C - rate?)

            #--------------------------------------------------------------------------------#
            # Concentration Calculations

            # Adding up where all of the substances are being consumed
            # diffusion + reaction + convection
            trypconc[i+1,:,:] = self.trypDiffusivity * (lapZArrays[0] + lapRArrays[0]) + tryp_rxn * crossSecAreas(rings, velocity_max)[0] + crossSecAreas(rings, velocity_max)[1] * (delZArrays[0] + delRArrays[0]) 
            serconc[i+1,:,:]  = self.serDiffusivity  * (lapZArrays[1] + lapRArrays[1]) + ser_rxn  * crossSecAreas(rings, velocity_max)[0] + crossSecAreas(rings, velocity_max)[1] * (delZArrays[1] + delRArrays[1]) 
            htpconc[i+1,:,:]  = self.htpDiffusivity  * (lapZArrays[2] + lapRArrays[2]) + htp_rxn  * crossSecAreas(rings, velocity_max)[0] + crossSecAreas(rings, velocity_max)[1] * (delZArrays[2] + delRArrays[2]) 



# Functions to calculate the Laplacian and deravities in r and z directions given our input matrix.
# Can be placed at the end of the script.
def laplacianZ(Z,dz):
    '''function to calculate second derivative by z'''
    Zleft = Z[1:-1,0:-2]
    Zright = Z[1:-1,2:]
    Zcenter = Z[1:-1,1:-1]
    return (Zleft - 2*Zcenter + Zright)/(dz**2)

def laplacianR(Z,r):
    '''Functions to calculate second derivative by r'''
    # Nate
    rArray = np.linspace(0,r,len(Z[0,:])) # what is this used for??
    dr = r/len(Z[0,:])
    Ztop = Z[0:-2,1:-1]
    Zbottom = Z[2:,1:-1]
    Zcenter = Z[1:-1,1:-1]
    return (Ztop - 2*Zcenter + Zbottom)/(dr**2) + deltaR(Z,r)/r

def deltaZ(Z,dz):
    '''function to calculate the first derivate by  z'''
    Zleft = Z[1:-1,0:-2]
    Zright = Z[1:-1,2:]
    return (Zright - Zleft)/(2*dz) 

def deltaR(Z,r):
    '''function to calculate the first derivate by  z'''
    #Nate and Marissa
    Ztop = Z[1:-1,0:-2]
    Zbottom = Z[1:-1,2:]
    dr = r/len(Z[0,:])
    return (Ztop - Zbottom)/(2*dr) 

def interpolateForValue(value, array):
    """
    Interpolation function designed to determine the M values for specific Graetz values

    Inputs:
    ------------------------------------------------------------------------------------
    value:
    array:

    Outputs:
    ------------------------------------------------------------------------------------
    interp: the interpreted variable
    """

    fromarray = array[:,0]
    for i in range(len(fromarray)):
        if value < fromarray[i]:
            break
    interp = []
    for n in range(10):
        toarray = array[:, n+1]
        interp.append(toarray[i]+((toarray[i]-toarray[i-1])/(fromarray[i]-fromarray[i-1]))*(value - fromarray[i]))
    return interp
                 
"""
    def getConcentration(self):
        
        A method designed to determine the bluk concentration in the intestines as it passes through the intestines.

        Outputs:
        ------------------------------------------------------------------------------------
        lengths: a list of length n that has the legnths at which the Graetz numbers were calculated at
        serGraetz: a list of Graetz numbers for ser (length n)
        trypGraetz: a list of Graetz numbers for tryp (length n)
        
        self.serMBetaList = serMBetaList = interpolateForValue(self.serEffPerm, self.MBeta)
        self.trypMBetaList = trypMBetaList = interpolateForValue(self.trypEffPerm, self.MBeta)
        for i in range(len(self.serConcentration)-1):

            # The initialize delta ser and tryp
            serDelta = 0
            trypDelta = 0

            # Calculation of delta ser and tryp over the first five terms of the series.
            for j in range(5):
                serDelta += serMBetaList[j+5]*np.exp(-1*serMBetaList[j]**2*self.serGraetz[i]) #Beta * 5 first in file, then M * 5
                trypDelta += trypMBetaList[j+5]*np.exp(-1*trypMBetaList[j]**2*self.trypGraetz[i])

            #REACTTTTTT
            c = (self.serConcentration[i], self.HTPConcentration[i], self.trypConcentration[i])
            serRate, HTPRate, trypRate = multistep(c, *self.kinetics)

            # Adding the concentration calculation to the list
            self.HTPConcentration[i+1] = self.HTPConcentration[i] + HTPRate*self.dt

            self.serConcentration[i+1] = self.serConcentration[i] * serDelta
            self.serConcentration[i+1] += serRate*self.dt

            self.trypConcentration[i+1] = self.trypConcentration[i] * trypDelta
            self.trypConcentration[i+1] += trypRate*self.dt
"""

