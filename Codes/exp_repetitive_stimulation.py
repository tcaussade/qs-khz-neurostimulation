from neuron import h
import numpy as np
from inh_src.functionalities import *

""""" 
Record action potentials for repetitive stimulation.
The folder to keep the files should be created before running this file.
"""""

def topulsearray(t,v, tr = -30):
    return pytools.detect_spikes(v, t, thresh=tr, fs= 20)

folder_for_saving = "insert_name_here"

# Search ranges
amplitudes = np.arange(0.0, 0.9, 0.015)
ftrains    = np.arange(0.1, 2.0+1e-5, 0.02)

# PARAMETERS
diameters     = np.array([7.3, 10.0, 12.8, 16.0])
fibermodel    = "MRG"
electricmodel = "IH"
active_duration = 40 # active duration of the source (ms)

# set waveform
forma = "canonical"
# forma = "asymmetric"
# forma = "shifted"

# set polarity
pol = +1
# pol = -1

# Select Media #
tissue_type = "brain_grey"
TISSUE = Tissue(tissue_type)
# Conductivity for quasi-static #
sqs = 0.105

# simulation resolution #
fmax = ftrains[-1] * 50
simulation_params = {"fcutoff" : fmax, # [kHz] 
                     "pp" : 25}        # Points to sample the shortest pulse
CompareExp = Experiment(simulation_params)
CompareExp.dt = 1e-6

for repetition_rate in ftrains:
    for amplitude in amplitudes:

        # source construction
        source_params = {"forma": forma, 
                        "amp": amplitude,
                        "ps": 10,  # ms
                        "pw": 0.2, # ms*2
                        "ftrain": repetition_rate, # kHz
                        "tend": 10 + active_duration + 10,
                        "npulses": int(active_duration * repetition_rate),
                        "loc": np.array([0.,1.0,0.]),
                        "pol": pol}

        # neuron parameters #
        h_params = {"temperature": 37, # Celsius
                    "v_init" : -77.3,  # [mV]
                    "tend" : source_params["tend"],
                    "cvode": 0}
        response = NeuronExperiment(h_params)

        # compute potentials & run neuron to compare #
        res_tot = dict()
        pots    = dict()

        source = ElectricPotential(source_params, electricmodel)
        t     = CompareExp.timegrid([source])
        pot   = CompareExp.compute_potential(t, source, TISSUE, sigma_qs = sqs)
        elec  = (source,pot)

        for fiberD in diameters:
            print("\nUsing fiberD %1.1f um" % fiberD)
            # create fiber #
            fiber_params = {"constant_cm" : 0,
                            "c_dc" : 1.0,
                            "is_xtra" : 1,
                            "intra_node" : 0,
                            "fiberD" : fiberD,
                            "nnodes" : 84}
            FIBER = FiberConstructor(fiber_params).assign_fiber_model(fibermodel)

            # instantiate neuron #
            fib_res = response.run_neuron(FIBER,t*1e3,(elec,))
            fib_spk = topulsearray(fib_res[0].as_numpy(),fib_res[1]["ve"])
            np.savetxt(folder_for_saving+"/spikestimestamp_ft%1.2f_amp%1.2f_fiberD%1.1f_%s.txt" % (source_params["ftrain"], amplitude,fiberD, electricmodel),
                    fib_spk, fmt = "%6g", delimiter=',', comments = '')
