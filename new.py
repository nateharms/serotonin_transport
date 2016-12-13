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

        self.serDiffusivity = serConditions.Diffusivity *3600 #6.2424e-8 m^2/sec
        self.serWallPerm = serConditions.Permeability*3600 #7.576e-13 m^2/sec
        self.trypDiffusivity = trypConditions.Diffusivity*3600 #5.386e-8 m^2/sec
        self.trypWallPerm = trypConditions.Permeability*3600 #6.44e-4 m^2/sec
        self.htpDiffusivity = htpConditions.Diffusivity*3600 #4.995e-8 m^2/sec
        self.htpWallPerm = htpConditions.Permeability *3600#estimating this as serotonin permeabiltiy - cannot find values for 5htp MLP


        self.permeabilities = (self.trypWallPerm,self.htpWallPerm,self.serWallPerm)
        self.diffusivities = (self.trypDiffusivity, self.htpDiffusivity, self.serDiffusivity)

        initalConcentrations = np.zeros(3, self.rings, self.sections)
        for i, J in enumerate(trypConditions, htpConditions, serConditions):
            initialConcentrations[i, :, 0] = J.Concentration
        times = np.arange(0, length/max_velocity, 10)
        result = odeint(self.getConcentration, initialConcentrations, times)

    def velocityProfile(self):
        r = np.linspace(0,self.radius,self.rings)
        velocityProfile = self.max_velocity*3600*(1-np.square(r)/np.square(self.radius))
        return velocityProfile

    def getConcentration(self, concentrationVector):
        rings = self.rings
        sections = self.sections
        velocityProf = self.velocityProfile()

        concentrations = concentrationsVector.reshape((3,N,M))
        trypConcentration = concentrationsVector[0]
        htpConcentration = concentrationsVector[1]
        serConcentration = concentrationsVector[2]

        concentrationList = [trypConcentration, htpConcentration, serConcentration]

        deltaZ = np.diff(concentrationsVector, axis = 1)
        deltaR = np.diff(concentrationsVector, axis = 2)
        lapZ = np.diff(concentrationsVector, n = 2, axis = 1)
        lapR = np.diff(concentrationsVector, n = 2, axis = 2)

        concentrationTuple = (trypConcentration, htpConcentration, serConcentration)

        #Rates from reaction kinetics
        tryp_rate, htp_rate, ser_rate = multistep(concentrationTuple, *self.kinetics)
        wallConcentrationTuple = (trypConcentration[ -1 , : ], htpConcentration[ -1 , : ], serConcentration[ -1 , : ])
        try_wallRate, htp_wallRate, ser_wallRate = multistep(wallConcentrationTuple, *self.wallKinetics)
        tryp_rate[ -1 , : ], htp_rate[ -1 , : ], ser_rate[ -1 , : ] = try_wallRate, htp_wallRate, ser_wallRate



        return rates.reshape(-1)
