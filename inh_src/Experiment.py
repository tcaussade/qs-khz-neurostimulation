import numpy as np
from alive_progress import alive_bar
from copy import deepcopy

##################################
# Instantiate single experiments #
##################################

class Experiment:

    def __init__ (self, params):
        self.fcutoff = params["fcutoff"] * 1e3 # Highest frequency harmonic
        self.pp      = params["pp"]            # Points to sample the shortest pulse
        self.dt      = np.Inf     
        self.skip_transient = True
        pass

    def display_ih_simulation_parameters(self, elec):
        self.ncoefs  = int(self.fcutoff * elec.tend) + 1
        print("\n Computing IH electric potential")
        print("Number of Fourier coefficients:", self.ncoefs - 1)
        print("Time step [ms]:", self.dt * 1e3)
        print("Number of points:", int(elec.tend/self.dt))

    def display_qs_simulation_parameters(self, sqs):
        print("\n Computing QS electric potential")
        print("Conductivity:", sqs, "[S/m]")
        
    # Time grid points #
    def timegrid(self, elec_list):
        tstart,tstop,tend = self.timegrid_params(elec_list)

        if self.skip_transient:
            d = 100*self.dt
            if tend > tstop:
                t_before = np.linspace(0.0, tstart - d, endpoint=False)
                t_during = np.arange(tstart-d,tstop+d, self.dt)
                t_after  = np.linspace(tstop+d+self.dt, tend)
            else:
                t_before = np.linspace(0.0, tstart - d, endpoint=False)
                t_during = np.arange(tstart-d,tend, self.dt)
                t_after  = np.array([])
        else:
            t_before = np.array([])
            t_during = np.arange(0.,tend, self.dt)
            t_after  = np.array([])

        return np.hstack((t_before, t_during, t_after)) 

    def timegrid_params(self, elecs_list):
        tstart  = np.Inf
        tstop   = 0.0
        tend    = 0.0
        for elec in elecs_list:
            tstart  = min(elec.ps, tstart)
            tstop   = max(elec.ps+elec.dur, tstop)
            tend    = max(elec.tend, tend)
            self.dt = min(elec.pw / self.pp, self.dt)
        return tstart,tstop,tend

    def _include_transient(self):
        self.skip_transient = False
    
    ##############################
    # Compute electric potential #
    ##############################

    def compute_ih_harmonic(self,elec,Tissue,fcoef,r,m):
        em,sm = Tissue.cole_cole_model(m*elec.wfun) 
        km    = m*elec.wfun * np.sqrt(Tissue.mu0*(em - 1j*sm/(m*elec.wfun)))
        sc    = ( sm + 1j*m*elec.wfun*em )
        phi   = elec.amp * fcoef(m) * np.exp(-1j*km * r) / (4*np.pi*sc * r)
        return phi 

    def compute_potential(self, t, elec, Tissue, r0 = np.array([0.,0.,0.]), sigma_qs = 0.0):
        
        r = np.sqrt( np.sum( (elec.loc - r0)**2 )  )

        if elec.model == "IH":
            self.display_ih_simulation_parameters(elec)
            fcoef, coefzero = elec.assign_fcoef()
            
            # m = 0 #
            phi   = elec.amp * coefzero() / (4*np.pi*Tissue.sigma_ionic * r) 
            tvec  = phi * np.ones(len(t), dtype = complex)

            # |m| > 0 #
            with alive_bar(self.ncoefs-1) as bar:
                for m in range(1,self.ncoefs):
                    bar()
                    phi   = self.compute_ih_harmonic(elec,Tissue,fcoef,r,m) * np.exp(1j*m*elec.wfun * t) 
                    tvec += phi + np.conj(phi)

        else:
            self.display_qs_simulation_parameters(sigma_qs)
            tvec = elec.qs_evaluate(t, elec.assign_qs() ) * elec.amp / (4*np.pi*sigma_qs * r)

        return elec.pol * np.real(tvec) * 1e3 # [mV]

    ##########################
    # Corrected conductivity #
    ##########################

    def optimal_sqs(self,tissue,elec_in):
        self.ncoefs = int(self.fcutoff * elec_in.tend) + 1
        elec = deepcopy(elec_in)
        elec.pw = 2 * elec_in.pw # DEBUG THIS LINE
        elec.n_shf = 1 # simplification
        fcoef, _ = elec.assign_fcoef()
        r = np.linalg.norm(elec.loc)

        # IH RMS-voltage
        eih = 0
        for m in range(1,self.ncoefs):
            em,sm = tissue.cole_cole_model(m*elec.wfun) 
            km    = m*elec.wfun * np.sqrt(tissue.mu0*(em - 1j*sm/(m*elec.wfun)))
            sc    = ( sm + 1j*m*elec.wfun*em )
            phim  = ( fcoef(m) * np.exp(-1j*km * r) / sc ) # / (4*np.pi*sc * r) 
            eih += 2 * np.abs(phim)**2 # |phim|**2 + |np.conj(phim)|**2
        eih *= 2*np.pi/elec.wfun

        # QS RMS-voltage
        eqs = elec.pw * elec.n_shf # / (4*np.pi*r)

        return np.sqrt(eqs/eih)