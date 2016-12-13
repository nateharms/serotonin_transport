import numpy as np
from matplotlib import pyplot as plt
from collections import namedtuple

import new

ConditionSet = namedtuple('ConditionSet', ['Concentration', 'Diffusivity', 'Permeability'])
KineticsParamers = namedtuple('ConditionSet', ['vmax1', 'Km1', 'K1', 'vmax2', 'Km2', 'K2'])


serCondition = ConditionSet(
                Concentration = 2,
                Diffusivity = 6.2424e-8,
                Permeability =7.576e-13)

trypCondition = ConditionSet(
                Concentration = .1,
                Diffusivity = 5.386e-8,
                Permeability = 6.44e-4)

htpCondition = ConditionSet(
                Concentration = 0,
                Diffusivity = 4.995e-8,
                Permeability = 7.576e-13)

kinetics = KineticsParamers(vmax1 = 0.1868167,
                            Km1 = 1.43533284,
                            K1 = 0.43680929,
                            vmax2 = 9.97704964,
                            Km2 = 2.37430847,
                            K2 = 0.25340153)

wallKinetics = KineticsParamers(vmax1 = 7.56e-14,
                            Km1 = 0.6e-3,
                            K1 = 1,
                            vmax2 = 7.14,
                            Km2 = 119,
                            K2 = 1)

radius = 2.5/2/100
length = 7.5
max_velocity = .0287/60
timestep = 20/3600
rings = 20
sections = 500




Model = new.LaminarFlow(length, radius, max_velocity, trypCondition, htpCondition, serCondition, kinetics, wallKinetics, rings, sections, timestep)

results = Model.results.reshape(-1,3,rings,sections)

np.save('results.npy', results)


# print(results.shape)

serotonin = results[:,2,:,:]
tryptophan = results[:,0,:,:]
htp = results[:,1,:,:]

serotoninUptake = np.array(Model.serotoninUptake)
np.save('serotoninUptake', serotoninUptake)
plt.plot(serotoninUptake)
plt.savefig('serotoninUptake.png', bbox_inches='tight')

from matplotlib import animation
import matplotlib
import seaborn as sns


matplotlib.rc('animation', html='html5')
fig = plt.figure()
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
timestamp = fig.text(0.45,0.1,'timestamp')
frames = 100
def animate(i):
    step = i * len(serotonin[:,0,0])//(frames)
    serotonin_g = serotonin[i,:,:50]
    tryp_g = tryptophan[i,:,:50]
    sns.heatmap(serotonin_g, vmin=0,vmax=1e-8, ax=ax1, cbar = None)
    ax1.set_title('Serotonin')
    sns.heatmap(tryp_g, vmin=0,vmax=1e-8, ax=ax2, cbar=None)
    ax2.set_title('Tryptophan')
    timestamp.set_text('t={0:.0f} hours'.format(Model.times[i]/3600))
    plt.tight_layout()



anim = animation.FuncAnimation(fig, animate, frames=frames, repeat_delay=2000, repeat=True)
anim
plt.savefig('dme_paper.png', bbox_inches='tight')
