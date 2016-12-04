import numpy as np

class LaminarFlow:
    def __init__(self, length, radius, serConcentration, trypConcentration, n, M_Beta_file = 'valuesBetaM.csv'):
        """
        Initializes the variables to run begin running the simulation: Maybe we make this general so we can use the same class
        to run Serotonin and Tryptophan. Haven't really decided yet. Different compounds would change our Graetz number as well.
        Beta and M values come from the mathematical solution to the differential equations.

        Input Variables
        ------------------------------------------------------------------------------------
        length: the length of the intestines (m)
        radius: the radius of the intestines (m)
        serConcentration: the initial concentration of seritonin in the intestines (mol/L)
        trypConcentration: the initial concentration of tryp in the intestines (mol/L)
        n: the number of sub-intervals of intestines looked at (unitless)
        M_Beta_file: the file that contains the M_Beta values needed for calculations
        """

        # Dimentions of the intestines
        self.length = length
        self.radius = radius

        # The initial concentrations for both of the tryp an ser.
        self.serConcentration = np.zeros(n)
        self.serConcentration[0] = serConcentration

        self.trypConcentration = np.zeros(n)
        self.trypConcentration[0] = trypConcentration

        # The diffusivities, wall permabilities, and effective permabilities for ser and tryp.
        self.serDiffusivity = 6.2424e-8 # m^2/sec
        self.serWallPerm = 7.576e-13 # m^2/sec ...?
        self.serEffPerm = self.serWallPerm * self.radius / self.serDiffusivity

        self.trypDiffusivity = 5.386e-8 # m^2/sec
        self.trypWallPerm = 6.44e-4 # m^2/sec
        self.trypEffPerm = self.trypWallPerm * self.radius / self.trypDiffusivity

        # The MBeta values
        self.MBeta = np.genfromtxt(M_Beta_file, delimiter = ',', skip_header = 1)

        self.getGraetz(n)
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
        lengths = []
        serGraetz = []
        trypGraetz = []
        max_velocity = 5
        for i in range(n):
            lengths.append(self.length * i / n)
            serGraetz.append(self.serDiffusivity*(lengths[i]/(max_velocity*self.radius**2)))
            trypGraetz.append(serGraetz[i]*self.trypDiffusivity/self.serDiffusivity)


        self.lengths = np.array(lengths)

        # The graetz numbers for both serotinin and trypophan.
        self.serGraetz = np.array(serGraetz)
        self.trypGraetz = np.array(trypGraetz)


    def getConcentration(self):
        self.hello = 'hello'
