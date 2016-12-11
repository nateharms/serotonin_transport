import numpy as np
from reactionKinetics import multistep


class LaminarFlow:
    '''
    Model for laminar flow through intestines, including reactions
    '''
    def __init__(self, length, radius, max_velocity, trypConditions, serConditions, htpConditions,
                    kinetics, wallKinetics, rings, sections, timestep): #M_Beta_file = 'valuesBetaM.csv'):
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
        htpConcentration: type: NamedTuple, with ('Concentration' (mM), 'Diffusivity', 'Wall Permeability')
        kinetics: (vmax1, Km1, K1, vmax2, Km2, K2) as variables for Michaelis Menten, hrs, mM in the bulk
        wallKinetics: (vmax1, Km1, K1, vmax2, Km2, K2) as variables for Michaelis Menten, hrs, mM at the wall
        rings: the number of ring partitions of intestines looked at (number of rings)
        sections: the number of sub-intervals of intestines looked at (number of sub-sections)
        timestep: duh (s)
        M_Beta_file: the file that contains the M_Beta values needed for calculations
        """

        # Dimentions of the intestines
        self.length = length
        self.radius = radius
        self.max_velocity = max_velocity
        self.kinetics = kinetics
        self.wallKinetics = wallKinetics

        residence_time = (length/max_velocity/3600)
        self.time = np.arange(0, residence_time, timestep)

        self.dt = timestep # to get hours

        # Subsection information
        self.rings = rings
        self.sections = sections

        # The initial concentrations for both of the tryp an ser.
        self.serConcentration = np.zeros((time, rings, sections))
        self.serConcentration[0,:,0] = serConditions.Concentration

        self.trypConcentration = np.zeros((time, rings, sections))
        self.trypConcentration[0,:,0] = trypConditions.Concentration

        self.htpConcentration = np.zeros((time, rings, sections))
        self.htpConcentration[0,:,0] = htpConditions.Concentration

        # The diffusivities, wall permabilities for tracked molecules
        # Effective permeability is no longer being used
        self.serDiffusivity = serConditions.Diffusivity #6.2424e-8 m^2/sec
        self.serWallPerm = serConditions.Permeability #7.576e-13 m^2/sec
        #self.serEffPerm = self.serWallPerm * self.radius / self.serDiffusivity

        self.trypDiffusivity = trypConditions.Diffusivity #5.386e-8 m^2/sec
        self.trypWallPerm = trypConditions.Permeability #6.44e-4 m^2/sec
        #self.trypEffPerm = self.trypWallPerm * self.radius / self.trypDiffusivity

        self.htpDiffusivity = htpConditions.Diffusivity #4.995e-8 m^2/sec
        self.htpWallPerm = htpConditions.Permeability #estimating this as serotonin permeabiltiy - cannot find values for 5htp MLP
        #self.htpEffPerm = self.serEffPerm

        self.permeabilities = (self.trypWallPerm,self.htpWallPerm,self.serWallPerm)
        self.diffusivities = (self.trypDiffusivity, self.htpDiffusivity, self.serDiffusivity)


        self.getConcentration()

    def getGraetz(self, n): #Not being used right now
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

    def crossSecAreas(self):
        """
        Determines the crossectional area of the rings based on given radius and number of subsections.

        Inputs:
        ------------------------------------------------------------------------------------
        self.rings: number of radial subsections chosen (unitless)
        self.radius: the radius of the intestines (m)

        Outputs:
        ------------------------------------------------------------------------------------
        totalArea: an array of areas corresponding to the number of rings chosen (m^2)
        velocityProfile: an array of velocities corresponding to radial position in the intestine (m/s)
        """

        totalArea = np.zeros(self.rings)
        velocityProfile = np.zeros(self.rings)
        r = self.radius/self.rings

        for m in range(rings):
            totarea[m] = np.pi * (r * (m+1)) ** 2 - np.pi * (r * m) ** 2
            velocityProfile[m] = velocity_max * (self.radius ** 2 - (r*m) **2 )

        return totalArea, velocityProfile

    def getConcentration(self):
        """
        This method is designed to determine the graetz number for both ser and tryp across the intestines.
        This value changes over the length of the intestines, however, not much else should.

        Inputs:
        ------------------------------------------------------------------------------------
        self.rings: number of radial subsections chosen (unitless)
        self.sections: number of subsections chosen in z direction (unitless)
        self.time: duh (sec)
        self.dt: the delta t for the system (sec?)
        self.trypConcentration: initial concentration profile of tryp (mM)
        self.serConcentration: initial concentration profile of ser (mM)
        self.htpConcentration: initial concentration profile of htp (mM)

        Outputs:
        ------------------------------------------------------------------------------------
        self.trypConcentration: final concentration profile of tryp (mM)
        self.serConcentration: final concentration profile of ser (mM)
        self.htpConcentration: final concentration profile of htp (mM)
        """
        # Creating three 3-D arrays to keeps track of the following:
        # ~ Substance Concentration
        # ~ Ring number (in the r direction)
        # ~ Section number (in the z direction)

        # Conventions that are used throughout the script, named for convienence.
        rings = self.rings
        sections = self.sections
        time = self.time # I am assuming this is an int

        # 3-D Arrays
        trypConcentration = self.trypConcentration
        serConcentration = self.serConcentration
        htpConcentration = self.htpConcentration

        ##### We pulled this from the self so Nate commented it out #####
        #trypConcentration = np.zeros((time, rings, sections))
        #serConcentration = np.zeros_like(trypConcentration)
        #htpConcentration = np.zeros_like(trypConcentration)

        # Adding the initial conditions to the tryp and ser arrays
        #trypConcentration[0,:,0] = trypConditions.concentration
        #serConcentration[0,:,0] = serConditions.concentration

        # Initializing the total amount of ser taken up.
        totalSerotoninUptake = 0

        # For loop to look through each time.
        for i in range(time-1):
            #--------------------------------------------------------------------------------#
            # Diffusion Calculations

            # Determining the laplacians and first deravitives for each concentration curve.
            # Laplicians are used for diffusion calculations
            # First deravitives are used for convection calculations

            delZArrays = []
            delRArrays = []
            lapZArrays = []
            lapRArrays = []

            concentrationTuple = (trypConcentration[i,:,:], serConcentration[i,:,:], htpConcentration[i,:,:])
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
            try_rxn = np.zeros_like(trypConcentration[i,:,:])
            ser_rxn = np.zeros_like(trypConcentration[i,:,:])
            htp_rxn = np.zeros_like(trypConcentration[i,:,:])

            # An iteration over the number of rings. I don't like this for loop, but it was the neatest way I though of how to deal with the reaction at the wall.
            try_rxn, ser_rxn, htp_rxn = reactionKinetics.multistep(concentrationTuple, *self.kinetics)

            #tryp_rxn, ser_rxn, htp_rxn = reactionKinetics.multistep(concentrationTuple, *self.kinetics) #trypConcentration, serConcentration, htpConcentration

            # Rxn at the surface
            #for Z in (tryp_rxn, ser_rxn, htp_rxn):
            #    Z[-1,:] *= 0 # Removing the existing reaction rates because they aren't correct for the surface.

            #for Z in concentrationTuple[0]:
            #    tryp_rxn[-1,:], ser_rxn[-1,:], htp_rxn[-1,:] = reactionKinetics.multistep(concentrationTuple)

            #--------------------------------------------------------------------------------#
            #Boundary Condition Layer

            #### I do not know what is going on here ####
            wallConcentrationTuple = (trypConcentrationZ[i,-1,:], serConcentrationZ[i,-1,:], htpConcentrationZ[i,-1,:])
            try_wallRxn, ser_rxn, htp_rxn = reactionKinetics.multistep(wallConcentrationTuple, *self.kinetics)

            wallDelta = [try_wallRxn, ser_rxn, htp_rxn]

            for j, Z in enumerate(trypConcentration, serConcentration, htpConcentration):
                Z[i,0,:] = Z[i,1,:]
                dcdr = self.permeabilities[j]/self.diffusivities[j](Z[i,-1,:]+wallDelta[j]*self.dt)
                Z[i,-1,:] = Z[i,-2,:]-dcdr*self.radius/rings

            #--------------------------------------------------------------------------------#
            # Concentration Calculations

            # Adding up where all of the substances are being consumed
            # diffusion + reaction + convection
            trypConcentration[i+1,:,:] = self.trypDiffusivity * (lapZArrays[0] + lapRArrays[0]) + tryp_rxn * crossSecAreas(rings, self.max_velocity)[0] + crossSecAreas(rings, self.max_velocity)[1] * (delZArrays[0] + delRArrays[0])
            serConcentration[i+1,:,:]  = self.serDiffusivity  * (lapZArrays[1] + lapRArrays[1]) + ser_rxn  * crossSecAreas(rings, self.max_velocity)[0] + crossSecAreas(rings, self.max_velocity)[1] * (delZArrays[1] + delRArrays[1])
            htpConcentration[i+1,:,:]  = self.htpDiffusivity  * (lapZArrays[2] + lapRArrays[2]) + htp_rxn  * crossSecAreas(rings, self.max_velocity)[0] + crossSecAreas(rings, self.max_velocity)[1] * (delZArrays[2] + delRArrays[2])

        # Assigning the final concentration arrays to self.
        self.trypConcentration = trypConcentration
        self.serConcentration = serConcentration
        self.htpConcentration = htpConcentration



# Functions to calculate the Laplacian and deravities in r and z directions given our input matrix.
# Can be placed at the end of the script.
def deltaZ(Z,dz):
    '''function to calculate the first derivative by z'''
    Zleft = Z[1:-1,0:-2]
    Zright = Z[1:-1,2:]
    return (Zright - Zleft)/(2*dz)

def deltaR(Z,r):
    '''function to calculate the first derivative by r'''
    Ztop = Z[1:-1,0:-2]
    Zbottom = Z[1:-1,2:]
    dr = r / len(Z[0,:])
    return (Ztop - Zbottom)/(2*dr)

def laplacianZ(Z,dz):
    '''function to calculate second derivative by z'''
    Zleft = Z[1:-1,0:-2]
    Zright = Z[1:-1,2:]
    Zcenter = Z[1:-1,1:-1]
    return (Zleft - 2*Zcenter + Zright)/(dz**2)

def laplacianR(Z,r):
    '''Functions to calculate second derivative by r'''
    dr = r / len(Z[0,:])
    Ztop = Z[0:-2,1:-1]
    Zbottom = Z[2:,1:-1]
    Zcenter = Z[1:-1,1:-1]
    first_der = deltaR(Z,r)
    rArray = np.linspace(0,r,len(first_der[:,0]))
    for i in range(len(first_der[:,0])):
        first_der[i,:] = rArray[i]*first_der[i,:]


    return (Ztop - 2*Zcenter + Zbottom)/(dr**2) + first_der

def interpolateForValue(value, array): # not being used right now
    """
    Interpolation function designed to determine the M values for specific Graetz values

    Inputs:
    ------------------------------------------------------------------------------------
    value: the value you are searching for
    array: the array you are searhing for 'value' in

    Outputs:
    ------------------------------------------------------------------------------------
    interp: the interpreted variable
    """

    fromarray = array[:,0]
    for i in range(len(fromarray)):
        if value < fromarray[i]:
            print('Could not find interpolated value. Out of range.')
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
