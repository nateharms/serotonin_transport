import numpy as np
from reactionKinetics import multistep

from time import sleep

class LaminarFlow:
    '''
    Model for laminar flow through intestines, including reactions
    '''
    def __init__(self, length, radius, max_velocity, trypConditions, htpConditions, serConditions,
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
        serConditions: type: NamedTuple, with ('Concentration' (mM), 'Diffusivity' (m^2/s), 'Wall Permeability'(m/s))
        trypConcentration: type: NamedTuple, with ('Concentration' (mM), 'Diffusivity' (m^2/s), 'Wall Permeability'(m/s))
        htpConcentration: type: NamedTuple, with ('Concentration' (mM), 'Diffusivity' (m^2/s), 'Wall Permeability'(m/s))
        kinetics: (vmax1, Km1, K1, vmax2, Km2, K2) as variables for Michaelis Menten, hrs, mM in the bulk
        wallKinetics: (vmax1, Km1, K1, vmax2, Km2, K2) as variables for Michaelis Menten, hrs, mM at the wall
        rings: the number of ring partitions of intestines looked at (number of rings)
        sections: the number of sub-intervals of intestines looked at (number of sub-sections)
        timestep: duh (hr)
        M_Beta_file: the file that contains the M_Beta values needed for calculations
        """

        #this is the main  variable we want to track
        # self.serotoninUptake = 0

        # Dimentions of the intestines
        self.length = length
        self.radius = radius
        self.max_velocity = max_velocity
        self.kinetics = kinetics
        self.wallKinetics = wallKinetics

        # Subsection information
        self.rings = rings
        self.sections = sections

        residence_time = (length/max_velocity/3600)
        self.time = np.arange(0, residence_time, timestep)
        self.serotoninUptake = np.zeros_like(self.time)
        self.dt = timestep # to get hours
        self.dz = length/sections


        self.outerRingSA = np.pi*2*radius*length/sections
        # The initial concentrations for both of the tryp an ser.
        self.serConcentration = np.zeros((len(self.time), rings, sections))
        self.serConcentration[0,:,0] = serConditions.Concentration

        self.trypConcentration = np.zeros((len(self.time), rings, sections))
        self.trypConcentration[0,:,0] = trypConditions.Concentration

        self.htpConcentration = np.zeros((len(self.time), rings, sections))
        self.htpConcentration[0,:,0] = htpConditions.Concentration

        # The diffusivities, wall permabilities for tracked molecules
        # Effective permeability is no longer being used
        self.serDiffusivity = serConditions.Diffusivity *3600 #6.2424e-8 m^2/sec
        self.serWallPerm = serConditions.Permeability*3600 #7.576e-13 m^2/sec
        #self.serEffPerm = self.serWallPerm * self.radius / self.serDiffusivity

        self.trypDiffusivity = trypConditions.Diffusivity*3600 #5.386e-8 m^2/sec
        self.trypWallPerm = trypConditions.Permeability*3600 #6.44e-4 m^2/sec
        #self.trypEffPerm = self.trypWallPerm * self.radius / self.trypDiffusivity

        self.htpDiffusivity = htpConditions.Diffusivity*3600 #4.995e-8 m^2/sec
        self.htpWallPerm = htpConditions.Permeability *3600#estimating this as serotonin permeabiltiy - cannot find values for 5htp MLP
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

    def velocityProfile(self):
        """
        Outputs:
        ------------------------------------------------------------------------------------
        velocityProfile: an array of velocities corresponding to radial position in the intestine (m/s)
        """

        r = np.linspace(0,self.radius,self.rings)
        velocityProfile = self.max_velocity*3600*(1-np.square(r)/np.square(self.radius))

        return velocityProfile

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

        rings = self.rings
        sections = self.sections
        time = self.time # I am assuming this is an int

        velocityProf = self.velocityProfile()

        # 3-D Arrays
        trypConcentration = self.trypConcentration
        serConcentration = self.serConcentration
        htpConcentration = self.htpConcentration

        for i in range(len(time)-1):
            # print('Running {} of {}'.format(i+1,len(time)-1))

            delZArrays = []
            delRArrays = []
            lapZArrays = []
            lapRArrays = []

            concentrationList = [trypConcentration[i,:,:], htpConcentration[i,:,:], serConcentration[i,:,:]]
            for j, Z in enumerate(concentrationList): #0 - tryp, 1 - htp, 2 - ser
                delZArrays.append(deltaZ(Z, self.dz))
                delRArrays.append(deltaR(Z, self.radius))
                lapZArrays.append(laplacianZ(Z, self.dz))
                # print('Length of laplacian is {}.'.format(len(lapZArrays[j])))
                lapRArrays.append(laplacianR(Z, self.radius))

            try_rxn = np.zeros_like(trypConcentration[i,:,:])
            ser_rxn = np.zeros_like(trypConcentration[i,:,:])
            htp_rxn = np.zeros_like(trypConcentration[i,:,:])

            concentrationTuple = (trypConcentration[i,:,:], htpConcentration[i,:,:], serConcentration[i,:,:])
            # An iteration over the number of rings. I don't like this for loop, but it was the neatest way I though of how to deal with the reaction at the wall.
            try_rxn, htp_rxn, ser_rxn = multistep(concentrationTuple, *self.kinetics)
            rxnList =  [try_rxn[1:-1,1:-1], htp_rxn[1:-1,1:-1], ser_rxn[1:-1,1:-1]]

            wallConcentrationTuple = (trypConcentration[i,-1,:], htpConcentration[i,-1,:], serConcentration[i,-1,:])
            try_wallRxn, htp_wallRxn, ser_wallRxn = multistep(wallConcentrationTuple, *self.wallKinetics)

            wallDelta = [try_wallRxn, htp_wallRxn, ser_wallRxn]



            #--------------------------------------------------------------------------------#
            # Concentration Calculations

            # Adding up where all of the substances are being consumed
            # diffusion + reaction + convection
            # + try_rxn[1:-2,1:-2]
            for j, Z in enumerate([trypConcentration, htpConcentration, serConcentration]):
                # print(len(delZArrays[j][:,0]))
                convection = np.zeros_like(Z[i+1,1:-1,1:-1])
                for k in range(len(velocityProf)-2):
                    # print('Length of first derivative is {} and length of velocity profile is {}.'.format(len(delZArrays[j][:,0]),len(velocityProf)))

                    assert len(delZArrays[j][:,0]) == len(velocityProf[1:-1])
                    # print(delZArrays[j])
<<<<<<< HEAD
                    convection[k,:] = velocityProf[k+1]*delZArrays[j][k,:]
=======
                    convection[k,:] = velocityProf[k]*delZArrays[j][k,:]
>>>>>>> 8d7229a790fb1d1d4f7700417f34e6dd6857b9a2
                diffusion = self.diffusivities[j]*(lapZArrays[j] + lapRArrays[j])
                Z[i+1,1:-1,1:-1] = Z[i,1:-1,1:-1] - ( convection )*self.dt #rxnList[j]



            # trypConcentration[i+1,1:-2,1:-2] = trypConcentration[i,1:-2,1:-2]+(self.trypDiffusivity * (lapZArrays[0] + lapRArrays[0])  + velocityProf[1:-2] * (delZArrays[0]))*self.dt
            # serConcentration[i+1,1:-2,1:-2]  = serConcentration[i,1:-2,1:-2]+(self.serDiffusivity  * (lapZArrays[2] + lapRArrays[2]) + ser_rxn[1:-2,1:-2]  + velocityProf[1:-2] * (delZArrays[2]))*self.dt
            # htpConcentration[i+1,1:-2,1:-2]  = htpConcentration[i,1:-2,1:-2]+(self.htpDiffusivity  * (lapZArrays[1] + lapRArrays[1]) + htp_rxn[1:-2,1:-2]  + velocityProf[1:-2] * (delZArrays[1]))*self.dt


            for j, Z in enumerate([trypConcentration, htpConcentration, serConcentration]):
                dc_thruwall = self.permeabilities[j]/self.radius*(Z[i,-1,:]) #+wallDelta[j]*self.dt
                #add diffusion term?
                Z[i+1,-1,1:-1] = Z[i,-1,1:-1] + (-1*dc_thruwall[1:-1]+self.diffusivities[j]*laplacianZ(Z[i,-1,:],self.dz+wallDelta[j][1:-1]))*self.dt
                if j == 2: #idk if this will work
                    # print(dcdr)
                    self.serotoninUptake[i+1] += self.serotoninUptake[i]+np.sum(dc_thruwall[1:-1]*self.radius/rings*self.outerRingSA)*self.dt
                #whatever the name of the boundary condition where the derivative is equal to 0
                Z[i+1,0,:] = Z[i+1,1,:]
                Z[i+1,:,0] = Z[i+1,:,1]
                Z[i+1,:,-1] = Z[i+1,:,-2]


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
    try:
        Zleft = Z[1:-1,0:-2]
        Zright = Z[1:-1,2:]
        Zcenter = Z[1:-1,1:-1]
    except IndexError:
        Zleft = Z[0:-2]
        Zright = Z[2:]
        Zcenter = Z[1:-1]
        # print('left:{}, right{}, center{}'.format(len(Zleft),len(Zright),len(Zcenter)))
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
        first_der[i,:] = first_der[i,:]/rArray[i]
    return (Ztop - 2*Zcenter + Zbottom)/(dr**2) + first_der
