import numpy as np

class LaminarFlow:
    def __init__(self, n, serConcentration, serGraetz, effPermeability, M_Beta_file = 'valuesBetaM.txt'):
        '''
        Initializes the variables to run begin running the simulation: Maybe we make this general so we can use the same class
        to run Serotonin and Tryptophan. Haven't really decided yet. Different compounds would change our Graetz number as well.
        Beta and M values come from the mathematical solution to the differential equations.
        '''
        self.graetz = serGraetz
        self.serConcentration = np.zeros(n)
        self.serConcentration[0] = serConcentration
        # self.trypConcentration = np.zeros(n)
        # self.trypConcentration[0] = trypConcentration
        self.effPermeability = effpermeability
        with open(M_Beta_file, 'r') as myfile:
            self.MBeta = np.fromfile(myfile)

        self.calculateConcentration()

    def calculateConcentration(self):
        MBetaList = interpolateForValue(self.effPermeability, self.MBeta)
        for i in range(len(serConcentration)):
            Delta = 0
            for j in range(5):
                Delta += MBetaList[j+5]*np.exp(MBetaList[j]**2*self.graetz) #Beta * 5 first in file, then M * 5

            serConcentration[i] = serConcentration[i-1]*Delta

def interpolateForValue(value, array):
    fromarray = array[0,:]
    for i in range(len(fromarray)):
        if value < fromarray[i]:
            break
    interp = []
    for n in range(10):
        toarray = array[n+1,:]
        interp.append(toarray[i]+((toarray[i]-toarray[i-1])/(fromarray[i]-fromarray[i-1]))*(value - fromarray[i]))
    return interp
