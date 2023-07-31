# from neuron import h
import numpy as np

from inh_src.functionalities import *

# PARAMETERS
ftrains         = np.arange(4,40.1,2)
electric_model  = "QS" # or "IH"
fibermodel      = "MRG"
diameters       = np.array([7.3, 10.0, 12.8, 16.0]) 
max_amp         = 1.5 # Maximal amplitude to perform binary-search 
active_duration = 80

wf       = "fixed" # or "fullduty", this refers to the block strategy
sqslabel = "sopt" # if corrected conductivity is desired
wflabel  = "_" + wf + "_" + sqslabel # file saving name

# Select Media #
tissue_type = "brain_grey"
TISSUE = Tissue(tissue_type)

# simulation resolution #
fmax = ftrains[-1] * 50
simulation_params = {"fcutoff" : int(fmax), # [kHz] ideally 500 kHz (Bossetti)
                    "pp" : 25}        # Points to sample the highest frequency
EXP = Experiment(simulation_params)
EXP.dt = 1e-7

for ft in ftrains:

    # Construct the source #
    source_params = {"forma": "biphasic_train", 
                    "amp": 1.0,
                    "ps": 10,  # ms
                    "pw": np.NaN,    # ms # 10kHz -> 0.1ms = 100us period
                    "ftrain": ft, # kHz
                    "tend": 10 + active_duration + 10,
                    "npulses": int(active_duration * ft),
                    "loc": np.array([0.,1.0,0.]),
                    "pol": -1}

    # pulse duration
    if wf == "fixed":
        source_params["pw"] = 0.01 * 2 # fixed pulsewidth at 10us
    elif wf == "fullduty":
        source_params["pw"] = 1/source_params["ftrain"] # full-duty cycle pulses

    # assemble the source
    elec  = ElectricPotential(source_params, electric_model)

    # quasi-static conductivity
    if sqslabel == "sopt":
        sqs  = EXP.optimal_sqs(TISSUE, elec)
    else: 
        sqs  = 0.105

    # compute potential
    t     = EXP.timegrid([elec])
    pot   = EXP.compute_potential(t,elec, TISSUE, sigma_qs = sqs)

    # neuron parameters #
    h_params = {"temperature": 37, # Celsius
                "v_init" : -77.3,  # [mV]
                "tend" : source_params["tend"],
                "cvode": 0}
    NeuronYale    = NeuronExperiment(h_params)

    # add intra_stimulation #
    intra_params = {"ps" : source_params["ps"]+40, "ftrain" : 1.0, "npulses" : 1}
    NeuronYale.set_intra_stim_on(intra_params)

    # Run neuron
    block_thresh = np.zeros(len(diameters))
    for i,fiberD in enumerate(diameters):
        # create fiber #
        fiber_params = {"constant_cm" : 0,
                        "c_dc" : 1.0,
                        "is_xtra" : 1,
                        "intra_node" : 0,
                        "fiberD" : fiberD,
                        "nnodes" : 84}
        FIBER = FiberConstructor(fiber_params).assign_fiber_model(fibermodel)   

        block_thresh[i] = NeuronYale.threshold_finder(FIBER, t*1e3, (elec,pot), amp1 = max_amp, tol= 1e-5)

        np.savetxt("inh_results/"+fibermodel+"/blocking"+wflabel+"/threshs_ft%1.2f_%s.txt" % (source_params["ftrain"], electric_model),
                        block_thresh, fmt='%6g', delimiter=',', comments='')
