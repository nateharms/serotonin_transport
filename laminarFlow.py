import numpy as np

class LaminarFlow:
    def __init__(self, length, radius, serConditions, trypConditions, n, M_Beta_file = 'valuesBetaM.csv'):
        """
        Initializes the variables to run begin running the simulation: Maybe we make this general so we can use the same class
        to run Serotonin and Tryptophan. Haven't really decided yet. Different compounds would change our Graetz number as well.
        Beta and M values come from the mathematical solution to the differential equations.

        Input Variables
        ------------------------------------------------------------------------------------
        length: the length of the intestines (m)
        radius: the radius of the intestines (m)
        serConditions: type: NamedTuple, with ('Concentration', 'Diffusivity', 'Wall Permeability')
        trypConcentration: type: NamedTuple, with ('Concentration', 'Diffusivity', 'Wall Permeability')
        n: the number of sub-intervals of intestines looked at (unitless)
        M_Beta_file: the file that contains the M_Beta values needed for calculations
        """

        # Dimentions of the intestines
        self.length = length
        self.radius = radius

        # The initial concentrations for both of the tryp an ser.
        self.serConcentration = np.zeros(n)
        self.serConcentration[0] = serConditions.Concentration

        self.trypConcentration = np.zeros(n)
        self.trypConcentration[0] = trypConditions.Concentration

        # The diffusivities, wall permabilities, and effective permabilities for ser and tryp.
        self.serDiffusivity = serConditions.Diffusivity #6.2424e-8 m^2/sec
        self.serWallPerm = serConditions.Permeability #7.576e-13 m^2/sec
        self.serEffPerm = self.serWallPerm * self.radius / self.serDiffusivity

        self.trypDiffusivity = trypConditions.Diffusivity #5.386e-8 m^2/sec
        self.trypWallPerm = trypConditions.Permeability #6.44e-4 m^2/sec
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

        max_velocity = 5
        self.lengths = lengths = np.linspace(0,self.length,n)
        self.serGraetz = self.serDiffusivity/(max_velocity*self.radius**2)*lengths
        self.trypGraetz = self.serGraetz*self.trypDiffusivity/self.serDiffusivity

    def getConcentration(self):
        """
        A method designed to determine the bluk concentration in the intestines as it passes through the intestines.

        Outputs:
        ------------------------------------------------------------------------------------
        lengths: a list of length n that has the legnths at which the Graetz numbers were calculated at
        serGraetz: a list of Graetz numbers for ser (length n)
        trypGraetz: a list of Graetz numbers for tryp (length n)
        """
        self.serMBetaList = serMBetaList = interpolateForValue(self.serEffPerm, self.MBeta)
        self.trypMBetaList = trypMBetaList = interpolateForValue(self.trypEffPerm, self.MBeta)
        for i in range(len(self.serConcentration)-1):

            # The initialize delta ser and tryp
            serDelta = 0
            trypDelta = 0

            # Calculation of delta ser and tryp over the first five terms of the series.
            for j in range(5):
                serDelta += serMBetaList[j+5]*np.exp(-1*serMBetaList[j]**2*self.serGraetz[i]) #Beta * 5 first in file, then M * 5
                trypDelta += trypMBetaList[j+5]*np.exp(trypMBetaList[j]**2*self.trypGraetz[i])

            # Adding the concetration calculation to the list
            self.serConcentration[i+1] = self.serConcentration[i] * serDelta
            self.trypConcentration[i+1] = self.trypConcentration[i] * trypDelta

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
