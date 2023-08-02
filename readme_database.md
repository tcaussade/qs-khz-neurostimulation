## Database

The folder `Database` contains all data utilized to make figures in 
*Towards a more accurate quasi-static approximation of the electric potential for neurostimulation with kilohertz-frequency sources* (T. Caussade, E. Paduro, M. Courdurier, E. Cerpa, W.E. Grill, L.E. Medina).


### Classical outputs.
Activation thresholds for fiber diameters of 7.3, 10, 12.8, and 16um. 
Different files employed for Current-Distance ("cd") and Strength-Duration ("sd"). 

Files are named as:
* `classical outputs/cd_thresholds_pwxxx_mod.txt`, where "xxx" is pulse duration (us).
* `classical outputs/sd_thresholds_distxxx_mod_yyy.txt`, where "xxx" is electrode-to-fiber distance (mm), and "mod" is the electric model (QS/IH). The conductivity choice "yyy" can be a fixed value ("fix") or the corrected conductivity ("sopt")
 

### Spike-train generated action potentials.
Timestamps of the generated action potentials using the "aaa" (canonical, asymmetric or shifted) biphasic waveform, or timestamps with canonical biphasic pulse train and anodic-first polarity has "aaa" as anodic

Files are named as `spiketimes_aaa/spiketimes_ftxxx_ampyyy_fiberDzzz_mod.txt`, where "xxx" is frequency (kHz), "yyy" is amplitude (mA), "zzz" is fiber diameter (um), and "mod" is the electric model (QS/IH).


### Conduction block thresholds

Block thresholds for fiber diameters of 7.3, 10, 12.8, and 16um, using fixed pulse durations of "aaa" um. If "aaa" is fullduty it refers to full-duty cycles, i.e. pulse duration is computed accordingly to the frequency of repetition.
The conductivity was "bbb", either using a fixed value ("fix") or the corrected conductivity ("corrected").

Files are named as `blocking_aaa_bbb/threshs_ftxxx_mod.txt`, where xxx is frequency (kHz), and "mod" is the electric model (QS/IH). 

Note the "corrected" folders do not contain any "IH" files, since the results are exactly the same as in "fix"
