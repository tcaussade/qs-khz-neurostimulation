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

# # Delete for outliers 
threshs_ih1[13,3] = np.NaN
threshs_ih2[16,2] = np.NaN

## "corrected" outliers
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

# plot full-duty cycle
fig, ax = plt.subplots(figsize = (12,8))

p1 = sns.lineplot(data = dfg, x = "freq", y = "fid",
                      hue = "fD", style = "model", palette = "flare", markers = True)
p1.set(ylim=(0,550), xlabel = "Frequency (kHz)", ylabel = "Block Threshold ($\mu$A)")
p1_err = p1.twinx()
sns.lineplot(data = errg, x = "freq", y = "err", ax = p1_err,
                      hue = "fD", style = "model", palette = "flare", markers = "^", linestyle = "dotted")
p1_err.set(ylim = (0,100), ylabel = "Percent error (%)")
h1,l1 = p1.get_legend_handles_labels()
h2,l2 = p1_err.get_legend_handles_labels()
l1.append("error"); h1.append(h2[-1])
l1 = ["7.3 $\mu$m", "10.0 $\mu$m", "12.8 $\mu$m", "16.0 $\mu$m", "QS", "IH", "percent error"]
h1.pop(5); h1.pop(0)

ax.legend(handles=h1, labels = l1, loc = 2)
p1_err.legend([],[], frameon = False)

# plot fixed duration
fig, ax2 = plt.subplots(figsize = (12,8))
p2 = sns.lineplot(data = dfs, x = "freq", y = "fid",
                      hue = "fD", style = "model", palette = "flare", markers = True)
p2.set(ylim = (0,1550), xlabel = "Frequency (kHz)", ylabel = "Block Threshold ($\mu$A)")
p2_err = p2.twinx()
sns.lineplot(data = errs, x = "freq", y = "err", ax = p2_err,
                      hue = "fD", style = "model", palette = "flare", markers = "^", linestyle = "dotted")
p2_err.set(ylim = (0,100), ylabel = "Percent error (%)")
h1,l1 = p1.get_legend_handles_labels()
h2,l2 = p1_err.get_legend_handles_labels()

l1.append("error"); h1.append(h2[-1])
l1 = ["7.3 $\mu$m", "10.0 $\mu$m", "12.8 $\mu$m", "16.0 $\mu$m", "QS", "IH", "percent error"]
h1.pop(5); h1.pop(0)

ax2.legend(handles=h1, labels = l1, loc = 1)
p2_err.legend([],[], frameon = False)

plt.show()