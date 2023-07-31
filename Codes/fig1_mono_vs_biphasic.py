# Monophasic and biphasic pulses getting closer to each other

import numpy as np
import matplotlib.pyplot as plt
from inh_src.ElectricPotential import *
from inh_src.Tissue import *
from inh_src.Experiment import *

tissue_type = "brain_grey"; TISSUE = Tissue(tissue_type)

source_params = {"forma": "biphasic_train", 
                "amp": 1.0,
                "ps": 1,  # ms
                "pw": 0.2,    # ms # 10kHz -> 0.1ms = 100us period
                "ftrain": .1, # kHz
                "tend": 200,
                "npulses": 5,
                "loc": np.array([0.,1.0,0.]),
                "pol": +1}

simulation_params = {"fcutoff" : 150.0, "pp" : 10} 
EXP = Experiment(simulation_params)
sqs = 0.105

active_dur = 10

ftrains = np.array([1/active_dur, .8, 2])

import seaborn as sns
import pandas as pd

potentials = dict()

for i,ft in enumerate(ftrains):
    for j,shp in enumerate(["monophasic_train", "biphasic_train"]):

        if shp == "biphasic_train": source_params["pw"] = 0.5
        else: source_params["pw"] = 0.25
        source_params["forma"]  = shp
        source_params["ftrain"] = ft
        source_params["npulses"]= int( active_dur * ft )

        elec_qs  = ElectricPotential(source_params, "QS")
        t_qs     = EXP.timegrid([elec_qs])
        pot_qs   = EXP.compute_potential(t_qs,elec_qs, TISSUE, sigma_qs = sqs)

        elec_ih  = ElectricPotential(source_params, "IH")
        t_ih     = EXP.timegrid([elec_ih])
        pot_ih   = EXP.compute_potential(t_ih,elec_ih, TISSUE)

        dfqs = pd.DataFrame({"t": t_qs*1e3, "pot":pot_qs, "model": "QS"})
        dfih = pd.DataFrame({"t": t_ih*1e3, "pot":pot_ih, "model": "IH"})

        potentials[i,j] = pd.concat([dfqs,dfih], keys=["QS","IH"])

tshow = 14

fig, ax = plt.subplots(2,len(ftrains), sharex = True, sharey = True, figsize = (12,6))
for i,ft in enumerate(ftrains):
    for j,shp in enumerate(["Monophasic", "Biphasic"]):
        data = potentials[i,j][potentials[i,j]["t"]<tshow]
        p = sns.lineplot(data = data, x = "t", y = "pot", hue = "model", ax = ax[j,i], palette=["#2aa5a5","#a52a2a"])
        p.set(xlabel = "Time (ms)", ylabel = shp+"\nAmplitude (mV)")
        ax[j,i].legend(loc = 4)
    ax[0,i].set_title( ("%1.0f Hz" % (ft*1e3)) )

        
plt.show()
