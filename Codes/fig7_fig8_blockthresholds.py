# Block thresholds

folder = "DataBase/blocking_"
sigma_lb = "corrected"
# sigma_lb = "fix"

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ftrains = np.arange(4,40.1,2)
diameters = ["7.3", "10.0", "12.8", "16.0"]

linestyles = {"QS" : "solid", "IH" : "dashed"}
colors     = {"7.3": "blue", "10.0" : "green", "12.8" : "orange", "16.0" : "red"}

# Extract data
threshs_qs1 = np.zeros( (len(ftrains), len(diameters)) )
threshs_ih1  = np.zeros( (len(ftrains), len(diameters)) )
threshs_qs2 = np.zeros( (len(ftrains), len(diameters)) )
threshs_ih2 = np.zeros( (len(ftrains), len(diameters)) )
for i,ft_data in enumerate(ftrains):
    ft = "%1.2f" % ft_data

    ######### Full duty cycle pulses ###########
    path_qs1 = folder+"fullduty_"+sigma_lb+ "/threshs_ft"+ft+"_QS.txt"
    threshs_qs1[i,:] = np.loadtxt(path_qs1, delimiter = ",")
    path_ih1 = folder+"fullduty_fix/threshs_ft"+ft+"_IH.txt"
    threshs_ih1[i,:] = np.loadtxt(path_ih1, delimiter = ",")

    ######### Fixed duration pulses #######
    path_qs2 = folder+"10_"+sigma_lb+ "/threshs_ft"+ft+"_QS.txt"
    threshs_qs2[i,:] = np.loadtxt(path_qs2, delimiter = ",")
    path_ih2 = folder+"10_fix/threshs_ft"+ft+"_IH.txt"
    threshs_ih2[i,:] = np.loadtxt(path_ih2, delimiter = ",")

# Data correction
threshs_ih1[13,3] = np.NaN
threshs_ih2[16,2] = np.NaN
if sigma_lb == "corrected":
    threshs_qs2[-4,1] = np.NaN
    threshs_qs2[12,0] = np.NaN
    threshs_qs2[16,2] = np.NaN

######### Errors ###########
errors1 = np.abs( (threshs_ih1-threshs_qs1)/threshs_ih1 ) * 100
errors2 = np.abs( (threshs_ih2-threshs_qs2)/threshs_ih2 ) * 100

def dataframe_creator(fid, model, info = "fid"):
    rep_ft   = np.repeat(ftrains, len(diameters))
    rep_diam = pd.to_numeric( np.array(diameters*len(ftrains)) )
    fid_size = len(ftrains)*len(diameters)        
    data = np.column_stack((rep_ft, rep_diam, fid.reshape(fid_size) * 1e3))
    cols = ["freq", "fD", info]
    df = pd.DataFrame(data, columns = cols)
    df["model"] = np.repeat(model, fid_size)
    return df

def assemble_experiment(fid_qs, fid_ih, exp):
    df_qs = dataframe_creator(fid_qs, "QS")
    df_ih = dataframe_creator(fid_ih, "IH")
    df    = pd.concat([df_qs, df_ih], keys = ["QS","IH"])
    df["exp"] = np.repeat(exp, 2*len(ftrains)*len(diameters))
    return df

dfg = assemble_experiment(threshs_qs1,threshs_ih1, "full-duty cycle")
dfs = assemble_experiment(threshs_qs2,threshs_ih2, "fixed pw")

def error_dataframe(df, exp):
    fid_qs = np.array(df[df["model"] == "QS"]["fid"])
    fid_ih = np.array(df[df["model"] == "IH"]["fid"])
    err    = np.abs(fid_qs - fid_ih) / np.abs(fid_ih) / 1e3 * 100
    return dataframe_creator(err, exp, info = "err")

errg = error_dataframe(dfg, "err")
errs = error_dataframe(dfs, "err")


plt.rc('font', size = 22)

err = []
for d in errg.freq.unique():
    err.append( errg.loc[errg.freq == d].err.max() )


# plot full-duty cycle
fig, ax = plt.subplots(figsize = (12,8))

p1 = sns.lineplot(data = dfg, x = "freq", y = "fid",
                      hue = "fD", style = "model", palette = "flare", markers = True, markersize = 10, dashes = [(1,0), (1,1)])
p1.set(ylim=(0,550), xlabel = "Frequency (kHz)", ylabel = r"Block Threshold ($\mu$A)")
p1_err = p1.twinx()
p1_err.plot(errg.freq.unique(), err, linestyle = "dashed", color = "black")
p1_err.set(ylim = (0,30), ylabel = "Percent error (%)")
h1,l1 = p1.get_legend_handles_labels()
h2,l2 = p1_err.get_legend_handles_labels()

l1  = [r"7.3 $\mu$m", r"10.0 $\mu$m", r"12.8 $\mu$m", r"16.0 $\mu$m"]
l1c = ["QS", "IH", "percent error"]
h1.pop(0)
h1c = []
h1c.append(h1[-1])
h1c.append(h1[-2])
h1c.append(h2[-1])
leg1 = ax.legend(handles=h1, labels = l1, loc = "upper left", bbox_to_anchor = (0, 1))
leg2 = ax.legend(handles=h1c, labels = l1c, loc = "upper left", bbox_to_anchor = (0.2, 1))
ax.add_artist(leg1)
p1_err.legend([],[], frameon = False)

err = []
for d in errs.freq.unique():
    err.append( errs.loc[errs.freq == d].err.max() )


# plot fixed duration
fig, ax2 = plt.subplots(figsize = (12,8))
p2 = sns.lineplot(data = dfs, x = "freq", y = "fid",
                      hue = "fD", style = "model", palette = "flare", markers = True, markersize = 10, dashes = [(1,0), (1,1)])
p2.set(ylim = (0,1550), xlabel = "Frequency (kHz)", ylabel = r"Block Threshold ($\mu$A)")
p2_err = p2.twinx()

p2_err.plot(errs.freq.unique(), err, linestyle = "dashed", color = "black")
p2_err.set(ylim = (0,30), ylabel = "Percent error (%)")
h1,l1 = p1.get_legend_handles_labels()
h2,l2 = p1_err.get_legend_handles_labels()

l1  = [r"7.3 $\mu$m", r"10.0 $\mu$m", r"12.8 $\mu$m", r"16.0 $\mu$m"]
l1c = ["QS", "IH", "percent error"]
h1.pop(0)
h1c = []
h1c.append(h1[-1])
h1c.append(h1[-2])
h1c.append(h2[-1])
leg1 = ax2.legend(handles=h1, labels = l1, loc = "upper right", bbox_to_anchor = (1, 1))
leg2 = ax2.legend(handles=h1c, labels = l1c, loc = "upper right", bbox_to_anchor = (0.8, 1))
ax2.add_artist(leg1)
p2_err.legend([],[], frameon = False)

plt.show()