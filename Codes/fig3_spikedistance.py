import numpy as np
import numpy as np
import matplotlib.pyplot as plt
import os

def spiketime_distance_spike(sa,sb,q): 
    # Spike train distance
    # Victor, Purpura. "Nature and precision of temporal coding in visual cortex: a metric-space analysis" (1996)
    m,n    = len(sa),len(sb)
    G      = np.zeros((m,n))
    G[:,0] = np.arange(m)
    G[0,:] = np.arange(n)
    for i in range(1,m):
        for j in range(1,n):
            G[i,j] = np.min([G[i-1,j]+1, G[i,j-1]+1, G[i-1,j-1]+q*np.abs(sa[i]-sb[j])])    
    return G[-1,-1]

folder = "DataBase/spiketimes_"
label  = "canonical"

amplitudes = np.arange(0.0, 0.9, 0.015)
ftrains    = np.arange(0.1, 2.0+1e-5, 0.02)
diameters  = ["7.3", "10.0", "12.8", "16.0"]

active_duration = 40 
stm_ps      = 10
pw          = 0.2

qvals = np.array([0,10,100])

spikedistance = np.zeros( (len(ftrains), len(amplitudes), len(qvals), len(diameters)) )

for m,fD in enumerate(diameters):
    for i, ft_data in enumerate(ftrains):
        ft = "%1.2f" % ft_data
        stm_nspikes = int(active_duration * ft_data)

        for j,amp_data in enumerate(amplitudes):
            amp = "%1.2f" % amp_data
            
            path_qs = folder+label+"/spikestimestamp_ft"+ft+"_amp"+str(amp)+"_fiberD"+str(fD)+"_QS.txt"
            path_ih = folder+label+"/spikestimestamp_ft"+ft+"_amp"+str(amp)+"_fiberD"+str(fD)+"_IH.txt"

            if os.path.exists(path_ih) and os.path.exists(path_qs):
                fib_spk_qs = np.loadtxt(path_qs); f1 = fib_spk_qs[fib_spk_qs>stm_ps]
                fib_spk_ih = np.loadtxt(path_ih); f2 = fib_spk_ih[fib_spk_ih>stm_ps]
                if (len(f1) == 0) or (len(f2)==0): 
                    continue
                else:
                    for k,q in enumerate(qvals):
                        spikedistance[i,j,k,m] = spiketime_distance_spike(f1, f2, q)

import matplotlib.colors as clr
cmap = clr.LinearSegmentedColormap.from_list('qs palette', ['#FFFFFF', "#a52a68"], N=256)

fig, ax = plt.subplots(len(qvals), len(diameters), sharex = True, sharey= True, figsize=(16,9))
X,Y = np.meshgrid(ftrains, amplitudes*1e3) 
for m,fD in enumerate(diameters):
    ax[0,m].set_title(fD+ " $\mu$m")
    ax[-1,m].set_xlabel("Frequency (kHz)")
    for k, q in enumerate(qvals):
        a = ax[k,m].scatter(X,Y, c = spikedistance[:,:,k,m].T, s = 30, cmap = cmap, vmax = 50)       
        fig.colorbar(a, ax = ax[k,m])
        ax[k,m].set_xticks([0,1/2,1,3/2,2])
for k,q in enumerate(qvals):
    ax[k,0].set_ylabel("q="+str(q)+ "\n Amplitude ($\mu$A)")

plt.rc('font', size = 18)
plt.show()
