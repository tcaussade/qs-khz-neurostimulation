import numpy as np
import matplotlib.pyplot as plt

import os

fibermodel = "MRG"
fD         = "12.8"

amplitudes = np.arange(0.0, 0.9, 0.015)
ftrains    = np.arange(0.1, 2.0+1e-5, 0.02)
assert len(ftrains)*len(amplitudes) <= 6000

active_duration = 40 # 80
stm_ps      = 10
pw          = 0.2


polarities = ["_biphasic_final", "_positive_polarity"]

fid_qs = np.zeros( (len(ftrains), len(amplitudes), len(polarities)) )
fid_ih = np.zeros( (len(ftrains), len(amplitudes), len(polarities)) )
err    = np.zeros( (len(ftrains), len(amplitudes), len(polarities)) )


for i, ft_data in enumerate(ftrains):
    ft = "%1.2f" % ft_data
    stm_nspikes = int(active_duration * ft_data)

    for j,amp_data in enumerate(amplitudes):
        for k, label in enumerate(polarities):
            amp = "%1.2f" % amp_data
            
            path_qs = "inh_results/"+fibermodel+"/spiketimes"+label+"/spikestimestamp_ft"+ft+"_amp"+str(amp)+"_fiberD"+str(fD)+"_QS.txt"
            path_ih = "inh_results/"+fibermodel+"/spiketimes"+label+"/spikestimestamp_ft"+ft+"_amp"+str(amp)+"_fiberD"+str(fD)+"_IH.txt"

            if os.path.exists(path_ih) and os.path.exists(path_qs):
                fib_spk_ih = np.loadtxt(path_ih)
                fib_spk_qs = np.loadtxt(path_qs)
            else:
                fib_spk_ih = np.array([])
                fib_spk_qs = np.array([])

            fid_qs[i,j,k] = len(fib_spk_qs[fib_spk_qs >= stm_ps])/stm_nspikes * 100
            fid_ih[i,j,k] = len(fib_spk_ih[fib_spk_ih >= stm_ps])/stm_nspikes * 100
            err[i,j,k]    = (fid_qs[i,j,k] - fid_ih[i,j,k]) / 100

# determine fcutoffs
fcut1 = np.zeros((len(amplitudes),len(polarities)))
fcut2 = np.zeros((len(amplitudes),len(polarities)))

fcut1_50 = np.zeros((len(amplitudes),len(polarities)))
fcut2_50 = np.zeros((len(amplitudes),len(polarities)))
for j in range(len(amplitudes)): # skip amp 0
    for k in range(len(polarities)):
        over_qs = ftrains[(fid_qs > 99.999)[:,j,k]]
        over_ih = ftrains[(fid_ih > 99.999)[:,j,k]]
        if np.any(over_qs):
            fcut1[j,k] = np.max( over_qs )
        if np.any(over_ih):
            fcut2[j,k] = np.max( over_ih )

import matplotlib.colors as clr
cmap_qs = clr.LinearSegmentedColormap.from_list('qs palette', ['#FFFFFF', '#2aa5a5'], N=256)
cmap_ih = clr.LinearSegmentedColormap.from_list('ih palette', ['#FFFFFF', '#a52a2a'], N=256)
cmap_err= clr.LinearSegmentedColormap.from_list('error palette', ['#a52a2a','#FFFFFF','#2aa5a5'], N=256)

X,Y = np.meshgrid(ftrains, amplitudes*1e3) 

fig, ax = plt.subplots(3,len(polarities), sharex = True, sharey= True, figsize = (8,9))
for k, lb in enumerate(polarities):
    a0 = ax[0,k].scatter(X,Y, c = fid_qs[:,:,k].T, s = 30, cmap = cmap_qs)
    a1 = ax[1,k].scatter(X,Y, c = fid_ih[:,:,k].T, s = 30, cmap = cmap_ih)
    if lb == "_positive_polarity":
        ax[0,k].set_title("Anodic")
    else: 
        ax[0,k].set_title("Cathodic")
    ax[-1,k].set_xlabel("Frequency (kHz)")

    cmax = np.max(np.abs(err)) # +0.05
    a2 = ax[2,k].scatter(X,Y, c = err[:,:,k].T, s = 30, cmap = cmap_err, vmin = -cmax, vmax = cmax)

    fig.colorbar(a0, ax = ax[0,k])
    fig.colorbar(a1, ax = ax[1,k])
    fig.colorbar(a2, ax = ax[2,k])

ax[0,0].set_ylabel("QS model \nAmplitude ($\mu$A)")
ax[1,0].set_ylabel("IH model \nAmplitude ($\mu$A)")
ax[2,0].set_ylabel("Errors \nAmplitude ($\mu$A)")
fig.suptitle("FiberD = "+fD)

plotcutoffs = True
if plotcutoffs:
    for k in range(len(polarities)):
        ax[0,k].plot(fcut1[:,k],amplitudes*1e3, color = "dimgrey", linestyle = "solid")
        ax[1,k].plot(fcut2[:,k],amplitudes*1e3, color = "black", linestyle = "dotted")
        ax[2,k].plot(fcut1[:,k],amplitudes*1e3, color = "dimgrey", linestyle = "solid")
        ax[2,k].plot(fcut2[:,k],amplitudes*1e3, color = "black", linestyle = "dotted")

# ax[0,0].set_xlim(0,2)
# ax[0,0].set_xticks(np.arange(0,2,0.5))
plt.show()
