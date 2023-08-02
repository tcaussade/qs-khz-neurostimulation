from neuron import h
import matplotlib.pyplot as plt
import numpy as np

from inh_src.functionalities import *

""""" 
Single strength-duration experiment to determine thresholds as function of pulse duration.
Experiment is repeated for the number of fiber diameters of interest.
"""""

# PARAMETERS
fiberdist  = 1.0  # electrode-to-fiber distance (mm)
model      = "QS" # electric model (QS/IH)
fibermodel = "MRG"
max_amp    = 0.5  # Maximal amplitude to perform binary-search 

diameters  = np.array([7.3, 10.0, 12.8, 16.0])
pulses_dur = 10 * np.logspace(-3,-1,20) 

sigma_qs   = "corrected" # determines the effective conductivity
# sigma_qs   = "naive"     # uses sqs = 0.105

# Select Media #
tissue_type = "brain_grey"
TISSUE = Tissue(tissue_type)

# simulation resolution #
simulation_params = {"fcutoff" : 500.0, # [kHz] 
                     "pp" : 30}         
StrengthDuration = Experiment(simulation_params)
StrengthDuration.dt = 1e-6

# neuron parameters #
h_params = {"temperature": 37, # Celsius
            "v_init" : -77.3,  # [mV]
            "tend" : 15,
            "cvode": 0}

# source # 
source_params = {"forma": "canonical", 
                "amp": 1.0,
                "ps": 10 ,     # [ms]
                "pw":  np.NaN, # [ms] 
                "ftrain": .1,  # [kHz]
                "tend": h_params["tend"],
                "npulses": 1,
                "loc": np.array([0.,fiberdist,0.]),
                "pol": -1}

diameters  = np.array([7.3, 10.0, 12.8, 16.0])
pulses_dur = 10 * np.logspace(-3,-1,20) 

thresh_current = np.zeros( (len(diameters), len(pulses_dur)), dtype = float)
for j,pw in enumerate(pulses_dur):
    source_params["pw"] = pw
    source = ElectricPotential(source_params, model)

    # Conductivity for quasi-static
    if sigma_qs == "corrected": # Correcetd value
        sqs = StrengthDuration.optimal_sqs(TISSUE, source) 
    elif sigma_qs == "naive": # fixed value
        sqs = 0.105 

    t      = StrengthDuration.timegrid([source])
    pot    = StrengthDuration.compute_potential(t, source, TISSUE, sigma_qs = sqs)

    for i,fiberD in enumerate(diameters):
        fiber_params = {"constant_cm" : 0,
                        "c_dc" : 1.0,
                        "is_xtra" : 1,
                        "intra_node" : 0,
                        "fiberD" : fiberD,
                        "nnodes" : 84}
        FIBER = FiberConstructor(fiber_params).assign_fiber_model(fibermodel)

        thresh_current[i,j] = NeuronExperiment(h_params).threshold_finder(FIBER, t*1e3, (source,pot), amp1  = max_amp)
