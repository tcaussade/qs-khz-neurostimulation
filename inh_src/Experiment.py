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

        # sampling_rate = self.pp * self.fcutoff
        # self.dt       = 1/sampling_rate
        self.dt       = np.Inf
        
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

    # def elec_timegrid(self, elec):
    #     return self.timegrid(elec.ps,elec.dur, elec.tend)

    # def timegrid(self, ps,dur,tend):
    #     # return np.arange(0.0, tend, self.dt)
    #     d = 20*self.dt
    #     if tend > ps+dur:
    #         t_before = np.linspace(0.0, ps - d, endpoint=False)
    #         t_during = ps + np.arange(-d,dur+d, self.dt)
    #         t_after  = np.linspace(ps+dur+d, tend)
    #     else:
    #         t_before = np.linspace(0.0, ps - d, endpoint=False)
    #         t_during = np.arange(ps-d,tend, self.dt)
    #         t_after  = np.array([])

    #     return np.hstack((t_before, t_during, t_after)) 

    def timegrid(self, elec_list):
        # return np.arange(0.0, tend, self.dt)
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
        # return t in [s]
        tstart  = np.Inf
        tstop   = 0.0
        tend    = 0.0
        for elec in elecs_list:
            tstart  = min(elec.ps, tstart)
            tstop   = max(elec.ps+elec.dur, tstop)
            tend    = max(elec.tend, tend)
            self.dt = min(elec.pw / self.pp, self.dt)
        # print("Changed dt [ms]:", self.dt * 1e3)
        return tstart,tstop,tend

    def _include_transient(self):
        self.skip_transient = False
    
    # Compute electric potential #

    def compute_ih_harmonic(self,elec,Tissue,fcoef,r,m):
        em,sm = Tissue.cole_cole_model(m*elec.wfun) 
        # sm    *= 1e-3 # sm should be in [S/mm]
        km    = m*elec.wfun * np.sqrt(Tissue.mu0*(em - 1j*sm/(m*elec.wfun)))
        sc    = ( sm + 1j*m*elec.wfun*em )
        phi   = elec.amp * fcoef(m) * np.exp(-1j*km * r) / (4*np.pi*sc * r)
        return phi 

    def compute_potential(self, t, elec, Tissue, r0 = np.array([0.,0.,0.]), sigma_qs = 0.0):
        
        # t = timegrid(elec.ps, elec.dur, elec.tend, dt = self.dt)
        # t = self.elec_timegrid(elec)
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
        # r is in [mm], amp is in [mA] and divide each other. 

    # # Fix Optimal QS conductivity

    def optimal_sqs(self,tissue,elec_in):
        self.ncoefs = int(self.fcutoff * elec_in.tend) + 1
        from copy import deepcopy
        elec = deepcopy(elec_in)
        elec.pw = 2 * elec_in.pw # WHY IS THIS ?
        elec.n_shf = 1 # safety: we set the integration for single pulses (assume the train has an unique pulsewidth)
        fcoef, _ = elec.assign_fcoef()
        r = np.linalg.norm(elec.loc)

        # IH energy
        eih = 0
        for m in range(1,self.ncoefs):
            em,sm = tissue.cole_cole_model(m*elec.wfun) 
            km    = m*elec.wfun * np.sqrt(tissue.mu0*(em - 1j*sm/(m*elec.wfun)))
            sc    = ( sm + 1j*m*elec.wfun*em )
            phim  = ( fcoef(m) * np.exp(-1j*km * r) / sc ) # / (4*np.pi*sc * r) 
            # eih +=  np.abs(phim)**2 + np.abs(np.conj(phim))**2 # |phim|**2 + |np.conj(phim)|**2
            # phim = self.compute_ih_harmonic(elec,tissue,fcoef,r,m)
            eih += 2 * np.abs(phim)**2
        eih *= 2*np.pi/elec.wfun

        # QS energy
        # eqs = elec.pw * elec.amp/ (4*np.pi*r)
        eqs = elec.pw * elec.n_shf

        return np.sqrt(eqs/eih)

    # def optimal_sqs(self,tissue,elec):
    #     assert elec.forma == "biphasic_train"
    #     etp    = deepcopy(elec)
    #     etp.ps = 1 * 1e-3
    #     etp.amp= 1.0 
    #     etp.tend = etp.ps + etp.pw * 2
    #     etp.n_shf= 1
    #     etp.pol  = +1
    #     etp.model = "IH"

    #     # source_params["pw"] = elec.pw * 1e3
    #     # source_params["loc"]= elec.loc * 1e3
    #     t     = self.timegrid([etp])
    #     pot   = self.compute_potential(t,etp, tissue) 
    #     pot   -=pot[0]  # Just for safety
    #     dist  = np.sqrt(np.sum(etp.loc**2))

    #     c_qs = 1.0 / (4*np.pi*dist) * (etp.pw) *1e3
    #     c_ih = self.quadrule( np.abs(pot) )  / dist * 1e-3 # np.sum(pot * EXP.dt) / dist
    #     return c_qs/c_ih 

    # def quadrule(self,y):
    #     # return self.dt * (y[0]/2 + np.sum(y[1:-2]) + y[-1]/2)
    #     return np.sum(y) * self.dt























# def polynomial_map(self,ta,tb):
#     n,p = [100,1]; c =  (n-1)**(p) 
#     t = np.linspace(ta,tb, endpoint=True, num = n)
#     tp = np.zeros(len(t))
#     for k in range(len(t)): tp[k] = ta + (tb-ta)/c * k**p
#     return tp

# def map_params(a,b,e1,e2):
# T = np.abs(b-a)
# c = 1 - np.log10(1-e2/T)/np.log10(e1/T)
# ki = 100 # k0
# err = 1
# i = 0
# maxiter = 15
# tol     = 1e-3
# f = lambda x: c*np.log10(x) - np.log10(x-1)
# df= lambda x: c/x - 1/(x-1) 
# # print("Newton Iterations:")
# while i < maxiter and err > tol:
#     ki = ki - f(ki)/df(ki)
#     err = np.abs(f(ki))
#     i +=1
#     # print(i,int(ki), np.round(err, decimals=5))
# p = np.ceil(- np.log10(e1/T)/np.log10(ki))
# return int(ki), int(p)

        # # t_before = self.polynomial_map(ps,0.0)[::-1]
        # t_before = np.linspace(0.0, ps - self.dt)
        # if ps+dur > tend: # useful for short simulations with single pulses.
        #     t_after = np.array([])
        #     # t_during = np.arange(ps + self.dt, tend, self.dt)
        #     t_during = ps + np.arange(-self.dt, tend, self.dt)
        # else: 
        #     # t_after  = self.polynomial_map(ps+dur,tend)
        #     # t_during = ps + np.arange(self.dt, dur, self.dt)
        #     t_during = ps + np.arange(-self.dt, dur+self.dt, self.dt)
        #     t_after  = np.linspace(ps+dur + self.dt, tend)

        # # print(np.abs(t_before[0]-t_before[1]), np.abs(t_before[-2]-t_before[-1]))
        # # print(t_during[0], t_during[-1])
        # # print(np.abs(t_after[0]-t_after[1]), np.abs(t_after[-2]-t_after[-1]))
        # return np.hstack((t_before, t_during, t_after))

