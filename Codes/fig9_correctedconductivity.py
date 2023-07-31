# Optimal value of sigma_qs as a function of pulse width

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from inh_src.Tissue import Tissue
from inh_src.Experiment import Experiment
from inh_src.ElectricPotential import ElectricPotential

tissues_type = ["brain_grey", "skin","fat", "muscle"]
waveforms = ["canonical", "asymmetric", "shifted"]

pws = np.logspace(-3,-0.4,20) * 2
df_sopt = pd.DataFrame()

simulation_params = {"fcutoff" : 2000.0, "pp" : 30}  
EXP = Experiment(simulation_params)
EXP.dt = 1e-6
EXP._include_transient()

list_dfs = []

for i, tissue in enumerate(tissues_type):
    TISSUE = Tissue(tissue)
    for j,wf in enumerate(waveforms):
        print("Computing for",tissue,"and",wf)

        ft = 1/(pws[-1]*2) / 2
        source_params = {"forma": wf, 
                "amp": 1.0, "ps": 0.5,"ftrain": ft, "tend": 10, "npulses": 1, "pol": -1, 
                "loc": np.array([0.,1.0,0.]), "pw": np.NaN}
        elec = ElectricPotential(source_params, "QS")

        ### Corrected SIGMA_QS ###
        sopt = np.zeros(len(pws))
        for n,pw in enumerate(pws):
            elec.pw = pw * 1e-3
            sopt[n] = EXP.optimal_sqs(TISSUE, elec)

        new_sopts = pd.DataFrame({"sopt" : sopt,
                                  "pw" : pws/2 * 1e3,
                                  "tissue": tissue,
                                  "wf" : wf}) 
        list_dfs.append(new_sopts)
df_sopt = pd.concat(list_dfs, keys = np.arange(len(list_dfs)) )


# plot
plt.figure(figsize=(4,3))
colors = ["#6b6867", "#875632"]
cmap = sns.color_palette("colorblind")
plt.rc('font', size = 16)

p = sns.lineplot(data = df_sopt, x = "pw", y = "sopt",
                 hue = "tissue", style = "wf", markers = True, palette = "colorblind")
ylims = np.arange(0.00,0.2501,0.01)
xlims = np.arange(0,250, 20)
p.set(xlabel = "Pulse duration ($\mu$s)", ylabel = "$\~\sigma_{qs}$ (S/m)",
      ylim = (ylims[0],ylims[-1]), yticks = ylims,
      xlim = (xlims[1]/10,xlims[-1]))
h,_ = p.get_legend_handles_labels()
l = ["Tissues", "Grey matter", "Skin", "Fat", "Muscle",
     "Waveforms","Canonical", "Asymmetric", "Shifted"]
p.legend(handles = h, labels = l)

from matplotlib.ticker import ScalarFormatter
ax = p.axes
ax.set_xscale("log")
ax.set_xticks([1,5,10,50,100])
ax.get_xaxis().set_major_formatter(ScalarFormatter())

plt.grid()
plt.show()