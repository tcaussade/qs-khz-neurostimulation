# Classic axonal outputs for biphasic pulses (Current-Distance and Strength-Duration)

folder = "DataBase/classical_outputs"

import numpy as np
import matplotlib.pyplot as plt

fibermodel = "MRG"
diameters  = ["7.3", "10.0", "12.8", "16.0"]

# pws       = ["50","100", "200", "500"]
pws = ["100"]
distances = 5 * np.logspace(-2,0,20)
cd_qs = np.zeros( (len(diameters), len(pws), len(distances)) )
cd_ih = np.zeros( (len(diameters), len(pws), len(distances)) )

pulses_dur = 5 * np.logspace(-3,-1,20) 
# distances   = ["0.05", "0.50", "1.00", "5.00"]
distances_sd = ["1.00"]
sd_qs = np.zeros( (len(diameters), len(distances_sd), len(pulses_dur)) )
sd_ih = np.zeros( (len(diameters), len(distances_sd), len(pulses_dur)) )

for fib_idx, fibD in enumerate(diameters):

    # Current-Distance data
    for idx,pw in enumerate(pws):
        cd_qs[fib_idx,idx, :] = np.loadtxt(folder+"/cd_thresholds_pw"+pw+"_"+"QS"+".txt", delimiter=",")[fib_idx,:]*1e3
        cd_ih[fib_idx,idx, :] = np.loadtxt(folder+"/cd_thresholds_pw"+pw+"_"+"IH"+".txt", delimiter=",")[fib_idx,:]*1e3

    # Strength-Duration data
    sd_explabel = "sopt" # "fix"
    for idx,dist in enumerate(distances_sd):
        sd_qs[fib_idx,idx, :] = np.loadtxt(folder+"/sd_thresholds_dist"+dist+"_"+"QS_"+sd_explabel+".txt", delimiter=",")[fib_idx,:]*1e3
        sd_ih[fib_idx,idx, :] = np.loadtxt(folder+"/sd_thresholds_dist"+dist+"_"+"IH_"+"fix"+".txt", delimiter=",")[fib_idx,:]*1e3

# Removing outliers (maximal binary-search ranges)
cd_qs[cd_qs == 1000] = np.NaN
cd_ih[cd_ih == 1000] = np.NaN
sd_qs[sd_qs == 500]  = np.NaN
sd_ih[sd_ih == 500]  = np.NaN

import pandas as pd
import seaborn as sns

# Building Current-Distance DF
rep_diam = pd.to_numeric( np.repeat(diameters, len(pws)*len(distances)) ) 
rep_pws  = pd.to_numeric( np.tile( np.repeat(pws, len(distances)), len(diameters) ))
rep_dists= np.tile(distances, len(pws)*len(diameters))
cd_size  = len(diameters)*len(pws)*len(distances)

data_cdqs = np.column_stack((rep_diam, rep_pws, rep_dists, cd_qs.reshape(cd_size)))
df_cdqs = pd.DataFrame(data_cdqs, columns=["fD", "pw", "dist", "thresh"])
df_cdqs["model"] = np.repeat("QS",cd_size)
data_cdih = np.column_stack((rep_diam, rep_pws, rep_dists, cd_ih.reshape(cd_size)))
df_cdih = pd.DataFrame(data_cdih, columns=["fD", "pw", "dist", "thresh"])
df_cdih["model"] = np.repeat("IH",cd_size)
df_cd = pd.concat([df_cdqs, df_cdih], keys = ["QS", "IH"])

err_cd = np.abs(cd_qs - cd_ih)/ np.abs(cd_ih) * 100
data_errcd = np.column_stack((rep_diam, rep_pws, rep_dists, err_cd.reshape(cd_size)))
df_errcd = pd.DataFrame(data_errcd,columns =  ["fD", "pw", "dist", "err"])
df_errcd["model"] = np.repeat("err",cd_size)

### Building Strength-Duration DF
rep_diam = pd.to_numeric( np.repeat(diameters, len(distances_sd)*len(pulses_dur)) ) 
rep_dists= pd.to_numeric( np.tile( np.repeat(distances_sd, len(pulses_dur)), len(diameters) ))
rep_pulse= np.tile(pulses_dur, len(distances_sd)*len(diameters)) * 1e3
sd_size  = len(diameters)*len(pulses_dur)*len(distances_sd)

data_sdqs = np.column_stack((rep_diam, rep_dists, rep_pulse, sd_qs.reshape(sd_size)))
df_sdqs = pd.DataFrame(data_sdqs, columns = ["fD", "dist", "pw", "thresh"])
df_sdqs["model"] = np.repeat("QS",sd_size)
data_sdih = np.column_stack((rep_diam, rep_dists, rep_pulse, sd_ih.reshape(sd_size)))
df_sdih = pd.DataFrame(data_sdih, columns = ["fD", "dist", "pw", "thresh"])
df_sdih["model"] = np.repeat("IH",sd_size)
df_sd = pd.concat([df_sdqs, df_sdih], keys = ["QS", "IH"])

err_sd = np.abs(sd_qs - sd_ih)/ np.abs(sd_ih) * 100
data_errsd = np.column_stack((rep_diam, rep_dists, rep_pulse, err_sd.reshape(sd_size)))
df_errsd = pd.DataFrame(data_errsd,columns =  ["fD", "dist", "pw", "err"])
df_errsd["model"] = np.repeat("err",sd_size)

import matplotlib
plt.rc('font', size = 22)

err = []
for d in df_errcd.dist.unique():
    err.append( df_errcd.loc[df_errcd.dist == d].err.max() )

# current-distance plot
ax_cd, ax = plt.subplots()
ax_cd = sns.lineplot(data = df_cd, x = "dist", y = "thresh", hue = "fD", style = "model", palette = "flare", markers = True, markersize = 10, dashes = [(1,0), (1,1)])
ax_cd.set(yscale = "log")
ax_errcd = ax_cd.twinx()
ax_errcd.plot(df_errcd.dist.unique(), err, linestyle = "dashed", color = "black")

ax_cd.set(xscale = "log", xlabel = "Distance (mm)", ylabel = r"Activation threshold ($\mu$A)") #, ylim=(0,500) )
ax_errcd.set(xscale = "log", ylim = (0,20), ylabel = "Percent error (%)")
h1,_ = ax_cd.get_legend_handles_labels()
h2,_ = ax_errcd.get_legend_handles_labels()

l1  = [r"7.3 $\mu$m", r"10.0 $\mu$m", r"12.8 $\mu$m", r"16.0 $\mu$m"]
l1c = ["QS", "IH", "percent error"]

h1.pop(0)
h1c = []
h1c.append(h1[-1])
h1c.append(h1[-2])
h1c.append(h2[-1])
leg1 = ax.legend(handles=h1, labels = l1, loc = "upper left", bbox_to_anchor = (0.0, 1))
leg2 = ax.legend(handles=h1c, labels = l1c, loc = "lower left", bbox_to_anchor = (0.0, 0.45))
ax.add_artist(leg1)
ax_errcd.legend([],[], frameon = False)

plt.show()

err = []
for d in df_errsd.pw.unique():
    err.append( df_errsd.loc[df_errsd.pw == d].err.max() )


# strength-duration plot
ax_sd, ax = plt.subplots()
ax_sd = sns.lineplot(data = df_sd, x = "pw", y = "thresh", 
                 hue = "fD", style = "model", palette = "flare", markers = True, markersize = 10, dashes = [(1,0), (1,1)])
ax_errsd = ax_sd.twinx()
ax_errsd.plot(df_errsd.pw.unique(), err, linestyle = "dashed", color = "black")
ax_sd.set(xscale = "log", xlabel = r"Pulse duration ($\mu$s)", ylabel = r"Activation threshold ($\mu$A)", ylim=(10,1000), yscale = "log")
ax_errsd.set(xscale = "log", ylim = (0,20), ylabel = "Percent error (%)")
h1,l1 = ax_sd.get_legend_handles_labels()
h2,l2 = ax_errsd.get_legend_handles_labels()

l1  = [r"7.3 $\mu$m", r"10.0 $\mu$m", r"12.8 $\mu$m", r"16.0 $\mu$m"]
l1c = ["QS", "IH", "percent error"]

h1.pop(0)
h1c = []
h1c.append(h1[-1])
h1c.append(h1[-2])
h1c.append(h2[-1])
leg1 = ax.legend(handles=h1, labels = l1, loc = "upper right", bbox_to_anchor = (1, 1))
leg2 = ax.legend(handles=h1c, labels = l1c, loc = "lower right", bbox_to_anchor = (1, 0.55))
ax.add_artist(leg1)
ax_errsd.legend([],[], frameon = False)

plt.show()


