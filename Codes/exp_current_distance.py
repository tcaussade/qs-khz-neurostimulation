from neuron import h
import numpy as np
from inh_src.functionalities import *

""""" 
Single current-distance experiment to determine thresholds as function of electrode-to-fiber distance.
Experiment is repeated for the number of fiber diameters of interest.
"""""

# PARAMETERS
pulse_width = 0.1 * 2 # pulse width (ms)
model       = "QS"    # electric model (QS/IH)
fibermodel  = "MRG"
max_amp     = 1.5     # Maximal amplitude to perform binary-search 

diameters = np.array([7.3, 10.0, 12.8, 16.0]) 
distances = 5 * np.logspace(-2,0,20)

# Select Media #
tissue_type = "brain_grey"
TISSUE = Tissue(tissue_type)

# Conductivity for quasi-static #
sqs = 0.105

# simulation resolution #
simulation_params = {"fcutoff" : 500.0, # [kHz] 
                     "pp" : 20}        # Points to sample the shortest pulse

CurrentDistance = Experiment(simulation_params)
CurrentDistance.dt = 1e-6

# neuron parameters #
h_params = {"temperature": 37, # Celsius
            "v_init" : -77.3,  # [mV]
            "tend" : 15,
            "cvode": 0}

# source # 
print("\nUsing %s model, pulse duration set to %1.0fus" % (model, pulse_width*1e3) )
source_params = {"forma": "biphasic_train", 
                "amp": 1.0,
                "ps": 10 ,     # [ms]
                "pw":  pulse_width, # [ms] 
                "ftrain": .1, # [kHz]
                "tend": h_params["tend"],
                "npulses": 1,
                "loc": np.array([0.,np.NaN,0.]),
                "pol": -1}

thresh_current = np.zeros( (len(diameters), len(distances)), dtype = float)

for i,fiberD in enumerate(diameters):

    fiber_params = {"constant_cm" : 0,
                    "c_dc" : 1.0,
                    "is_xtra" : 1,
                    "intra_node" : 0,
                    "fiberD" : fiberD,
                    "nnodes" : 84}
    FIBER = FiberConstructor(fiber_params).assign_fiber_model(fibermodel)

    for j,dist_to_fiber in enumerate(distances):
        source_params["loc"] = np.array([0.,dist_to_fiber,0.])
        source = ElectricPotential(source_params, model)
        t      = CurrentDistance.timegrid([source])
        pot    = CurrentDistance.compute_potential(t, source, TISSUE, sigma_qs = sqs)

        thresh_current[i,j] = NeuronExperiment(h_params).threshold_finder(FIBER, t*1e3, (source,pot), amp1 = max_amp)
