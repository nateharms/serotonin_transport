class LaminarFlow:
    def __init__(self, serConcentration, trypConcentration, graetz, effPermeability):
        self.graetz = graetz
        self.serConcentration = serConcentration
        self.trypConcentration = trypConcentration
        self.effPermeability = effpermeability

    def calculateConcentration(self):
        #insert differential equations

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
