import numpy as np
import matplotlib.pyplot as plt
import os

fibermodel = "MRG"
label      = "_biphasic_final"

amplitudes = np.arange(0.0, 0.9, 0.015)
ftrains    = np.arange(0.1, 2.0+1e-5, 0.02)
diameters  = ["7.3", "10.0", "12.8", "16.0"]

active_duration = 40 # 80
stm_ps      = 10
pw          = 0.2

# Fix values
idx_amp = np.searchsorted(amplitudes, 0.405)
amp_ref = "%1.2f" % amplitudes[idx_amp]
idx_ft  = np.searchsorted(ftrains, [0.3, 0.76, 1.2, 1.9])

ft_lb = []
for idx in idx_ft:
    if ftrains[idx] < 1:
        ft_ref  = "%1.0f Hz" % (ftrains[idx]*1e3)
    else:
        ft_ref  = "%1.1f kHz" % (ftrains[idx])
    ft_lb.append(ft_ref)

idx_fib = np.searchsorted(diameters, "12.8")
fD      = diameters[idx_fib]

fig, ax = plt.subplots(len(idx_ft),1, sharex=True, figsize = (8,4))
i = 0
for j,idx in enumerate(idx_ft):
    ft_ref  = "%1.2f" % ftrains[idx]
    # Read the data
    path_qs = "inh_results/"+fibermodel+"/spiketimes"+label+"/spikestimestamp_ft"+ft_ref+"_amp"+str(amp_ref)+"_fiberD"+str(fD)+"_QS.txt"
    path_ih = "inh_results/"+fibermodel+"/spiketimes"+label+"/spikestimestamp_ft"+ft_ref+"_amp"+str(amp_ref)+"_fiberD"+str(fD)+"_IH.txt"
    if os.path.exists(path_ih) and os.path.exists(path_qs):
        fib_spk_qs = np.loadtxt(path_qs); f1 = fib_spk_qs[fib_spk_qs>stm_ps]
        fib_spk_ih = np.loadtxt(path_ih); f2 = fib_spk_ih[fib_spk_ih>stm_ps]
    # source stamps
    nsource = int(active_duration * float(ft_ref))
    source = stm_ps + np.arange(nsource)/float(ft_ref)

    # Raster plot
    ax[j].vlines(f1, 1+i, 2+i, color = "#2aa5a5", label = "QS", lw = 3)
    ax[j].vlines(f2, 0+i, 1+i, color = "#a52a2a", label = "IH", lw = 3)
    ax[j].vlines(source, 0.75+i,1.25+i, color = "grey", label = "source")

    # i+= 3

    ax[j].set_xlim([stm_ps-5, stm_ps+active_duration+5])
    ax[j].set_yticks([])
    ax[j].set_ylabel(ft_lb[j])

ax[0].legend(["QS", "IH", "Source"], loc = 0)
ax[-1].set_xlabel('Time (ms)')
plt.show()