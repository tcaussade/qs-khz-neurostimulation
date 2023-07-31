from neuron import h
import numpy as np

import inh_src.pytools.tools as pytools

class FiberConstructor:

    def __init__(self, params):

        self.constant_cm = params["constant_cm"] # 0: freq-dependent cm
        self.c_dc = params["c_dc"]
        self.fiberD = params["fiberD"]
        self.nnodes = params["nnodes"]
        self.is_xtra= params["is_xtra"]
        pass

    def assign_fiber_model(self, model):
        h.load_file("nrnfiberlib/fiber_library.hoc")
        h.load_file('stdrun.hoc')
        self.fibermodel = model
        if model == "MRG":
            return  h.MRGFiber(self.fiberD, self.nnodes)
        elif model == "DCFiber":
            return  h.DCFiber(self.fiberD, self.nnodes,
                       self.constant_cm, self.c_dc, self.is_xtra)
        pass

    def assign_neuron_model(self, model):
        pass

class NeuronExperiment:

    def __init__(self, params):
        self.temperature = params["temperature"] # 37    # Celsius degrees
        self.v_init      = params["v_init"]      # -77.3 # [mV]
        self.tend        = params["tend"]
        self.cvode       = params["cvode"]

        self.intra_stim = False
        self.recordall  = False
        self.block_test = False
        pass

    def set_intra_stim_on(self, intra_params):
        self.intra_stim = True
        self.intra_ps   = intra_params["ps"]
        self.intra_ft   = intra_params["ftrain"]
        self.intra_nshf = intra_params["npulses"]
        self.block_test = True
        print("threshold_finder() method now computes block thresholds")
        pass
    def set_intra_stim_off(self):
        self.intra_stim = False
        pass

    def record_allnodes_on(self):
        self.recordall = True
        pass
    def record_allnodes_off(self):
        self.recordall = False
        pass


    def run_neuron(self,fiber,te,pot_list, verbose = True):

        tstim = h.Vector(len(te))
        tstim.from_python(te)
        dt_neuron = np.min(np.diff(te)) # Neuron simulation time_step

        # Set extracellular stimuls
        extra = dict()
        if verbose: print("Applying extracellular potential...")
        for i,s in enumerate(fiber.sl):
            locf = np.array([fiber.xcoord[i],fiber.ycoord[i],fiber.zcoord[i]])
            local_pot = np.zeros(len(te))
            for elec,pot in pot_list:
                rf        = np.sqrt( np.sum((locf - elec.loc*1e3)**2) ) 
                r0        = np.sqrt( np.sum((elec.loc*1e3)**2) ) # scale to fiber distance 
                local_pot += np.real(pot) * r0 / rf # I/r0 * r0/r.loc = pot * r0/r.loc
                # if i == midnode: print("dist is %.4f" % rf) # uncomment to check scaling factor by distance
            extra[i] = h.Vector()
            extra[i].from_python(local_pot) # stimulation from all electrodes
            extra[i].play(s(.5)._ref_e_extracellular,tstim, True)

        # Set intracellular stimulus
        if self.intra_stim: # used to test signal blocking
            if verbose: print("Adding intracellular stimulus...")
            intra_t, intra_svec = self.intra_stimulus()  # Intra-stimulation node
            iclamp = h.IClamp(.5, sec=fiber.node[1])
            iclamp.delay = 0
            iclamp.dur = 1e9
            iclamp.amp = 0.1
            intra_svec.play(iclamp._ref_amp, intra_t, 1)
 

        h.celsius = self.temperature
        h.v_init  = self.v_init
        h.tstop   = self.tend
        h.dt = dt_neuron

        cvode = h.CVode(); cvode.active(0)         # CVode tools
        #cvode.use_daspk(1)
        #h.secondorder = 2 # NEURON: NrnDAEs only work with secondorder==0 or daspk

        # Outputs #
        res = dict()

        endnode = int(fiber.nnodes - 2) 
        ve = h.Vector(); ve.record(fiber.node[endnode](.5)._ref_v, tstim)
        res["ve"] = ve # Potential at the end of the fiber # 

        if self.recordall: # save every node (if True)
            vnodes = dict() 
            snodes = dict()
            for i in range( int(fiber.nnodes) ):
                v = h.Vector()
                v.record(fiber.node[i](.5)._ref_v, tstim)
                vnodes[i] = v

                s = h.Vector()
                s.record(fiber.node[i](.5)._ref_e_extracellular, tstim)
                snodes[i] = s
            res["allnodes"] = vnodes
            res["allextra"] = snodes

        midnode = int(fiber.nnodes / 2) 
        ex = h.Vector(); ex.record(fiber.node[midnode](.5)._ref_e_extracellular, tstim)
        res["ex"] = ex # Extracellular potential at the center of the fiber #
        
        th = h.Vector(); th.record(h._ref_t, tstim) # time #

        
        if verbose: print("Running neuron!")
        h.init() # Run simulation!
        h.run() 
        
        return th, res

    def ask_activation(self,t,ve):
        spk = pytools.detect_spikes(ve, t.as_numpy(), thresh=-30)
        return np.any(spk)
    def ask_blocking(self,t,ve):
        endfiber = pytools.detect_spikes(ve, t.as_numpy(), thresh=-30)
        return not (np.any(endfiber[endfiber<self.intra_ps+3]>self.intra_ps))

        # assert len(endfiber) < 3 # check there are no more than 1 pulse (without considering the first one)
        return not ( np.any(endfiber>self.intra_ps) * np.any(endfiber<self.intra_ps+3) )
    def question(self,t,ve , to = 0.0):
        if self.block_test:
            return self.ask_blocking(t,ve)
        else: # activation_test
            return self.ask_activation(t,ve)

    def threshold_finder(self,fiber,t, source, amp0 = 0.0, amp1 = 1.0, tol = 1e-8):
        source_tp = np.copy(np.asarray(source, dtype = object))

        iter_count = 0
        while amp1-amp0 > tol:
            amp = 0.5 * (amp1+amp0)
            source_tp[1] = amp * source[1] # source = [info,pot]

            th,res = self.run_neuron(fiber, t, (source_tp,), verbose = False)
            # print("testing for I =", amp)
            
            if self.question(th,res["ve"]): 
                amp1 = amp
            else: 
                amp0 = amp
            iter_count +=1
        print("Converged in %i iterations" % iter_count )
        print("\nThreshold found at: \nA = %.4f [mA]\n" % amp)
        return amp
    
    def intra_stimulus(self):

        from inh_src.ElectricPotential import ElectricPotential
        from inh_src.Experiment import Experiment

        intra_params = {"forma": "monophasic_train", 
                        "amp": 1.0,
                        "ps": self.intra_ps,  # ms
                        "pw": 0.03,    # ms # 10kHz -> 0.1ms = 100us period
                        "ftrain": self.intra_ft , # kHz
                        "tend": self.tend,
                        "npulses": self.intra_nshf,
                        "loc": np.array([0.,1.,0.]),
                        "pol": +1}

        intra_sim_params =  {"fcutoff" : 50.0, "pp" : 10} 
        intra_elec = ElectricPotential(intra_params, "QS")
        Experiment(intra_sim_params)._include_transient()
        intra_t    = Experiment(intra_sim_params).timegrid([intra_elec])
        intra_pot  = Experiment(intra_sim_params).compute_potential(intra_t,intra_elec, np.NaN, sigma_qs = 0.105)
        
        intrastimvec = h.Vector(len(intra_pot))
        intrastimvec.from_python(intra_pot)

        # import matplotlib.pyplot as plt
        # plt.plot(intra_t, intra_pot)
        # plt.plot(intra_t, intrastimvec)
        # plt.title("intra stimulus")
        # plt.show()

        intra_tstim = h.Vector(len(intra_t))
        intra_tstim.from_python(intra_t * 1e3)


        return intra_tstim, intrastimvec