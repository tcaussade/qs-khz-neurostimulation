import numpy as np
from numpy import pi

from alive_progress import alive_bar

class ElectricPotential:

    def __init__(self, params, model):
        # Units are in SI
        self.forma = params["forma"]
        self.ps    = params["ps"] * 1e-3 # [s]
        self.pw    = params["pw"] * 1e-3 # [s]
        self.tend  = params["tend"] * 1e-3 # [s]
        self.ftrain= params["ftrain"] * 1e3 # [Hz]
        self.n_shf = params["npulses"]
        self.amp   = params["amp"] * 1e-3 # [A]

        self.loc   = params["loc"] * 1e-3 # [m]
        self.pol   = params["pol"]

        self.ffreq  = 1/self.tend
        self.dur    = 1/self.ftrain * self.n_shf

        self.wfun  = 2*np.pi/self.tend 
        self.wtrain= 2*np.pi * params["ftrain"]

        self.model = model # QS or IH
        pass

    #########################
    # Fourier coefficients
    #########################

    def assign_fcoef(self):
        if self.forma == "monophasic":
            return self.fcoef_monophasic, self.coefzero_monophasic
        elif self.forma == "canonical":
            return self.fcoef_canonical, self.coefzero_energy
        elif self.forma == "asymmetric":
            return self.fcoef_asymmetric, self.coefzero_energy
        elif self.forma == "shifted":
            return self.fcoef_shifted, self.coefzero_energy
        else:
            print("NO FCOEF HAS BEEN SELECTED")
            print("\ncheck your spelling...")

    def coefzero_energy(self,*args):
        # Biphasic pulses in general
        return 0.0

    ####### Train of *n_shf* monophasic pulses
    def fcoef_monophasic(self, m):
        c = 0
        for n in range(self.n_shf):
            p1 = self.ps + n/self.ftrain
            p2 = p1 + self.pw
            c += np.exp(-1j*m*self.wfun*p2) - np.exp(-1j*m*self.wfun*p1)
        return c * 1j/(2*pi*m)

    def coefzero_monophasic(self):
        return self.n_shf * self.pw * self.ffreq

    ####### Train of *n_shf* biphasic pulses
    def fcoef_canonical(self, m):
        c = 0
        for n in range(self.n_shf):
            p1 = self.ps + n/self.ftrain
            p2 = p1 + self.pw
            pm = 0.5 * (p1+p2)
            c += - np.exp(-1j*m*self.wfun*p2) + 2*np.exp(-1j*m*self.wfun*pm) - np.exp(-1j*m*self.wfun*p1)
        return c * 1j/(2*pi*m)
    
    def fcoef_asymmetric(self,m):
        c = 0
        for n in range(self.n_shf):
            p1 = self.ps + n/self.ftrain
            pm = p1 + self.pw * 0.5
            pe = p1 + self.pw * 2.5
            c += 1.25*np.exp(-1j*m*self.wfun * pm) - np.exp(-1j*m*self.wfun * p1) - 0.25*np.exp(-1j*m*self.wfun * pe)
        return c * 1j/(2*np.pi*m)
    
    def fcoef_shifted(self,m):
        c = 0
        for n in range(self.n_shf):
            p1  = self.ps + n/self.ftrain
            pe1 = p1 + self.pw * 0.5
            pe2 = p1 + 0.5/self.ftrain
            pe3 = pe2 + 0.5*self.pw
            c += np.exp(-1j*m*self.wfun * pe1) - np.exp(-1j*m*self.wfun * p1) - np.exp(-1j*m*self.wfun * pe3) + np.exp(-1j*m*self.wfun * pe2)
        return c * 1j/(2*np.pi*m)
 
    #########################
    # Quasi-static model
    #########################

    def assign_qs(self):
        if self.forma == "monophasic":
            return self.qs_monophasic
        elif self.forma == "canonical":
            return self.qs_canonical 
        elif self.forma == "asymmetric":
            return self.qs_asymmetric
        elif self.forma == "shifted":
            return self.qs_shifted

    def qs_monophasic(self, t, n):
        p1 = self.ps + n/self.ftrain
        p2 = p1 + self.pw
        return np.piecewise(t,[t<=p1,(p1<t)*(t<p2),t>=p2],[0.0,1.0,0.0])  

    def qs_canonical(self, t, n):
        p1 = self.ps + n/self.ftrain
        p2 = p1 + self.pw
        pm = 0.5 * (p1+p2)
        return np.piecewise(t,[t<=p1,(p1<t)*(t<pm),t>=pm],[0.0,+1.0,0.0]) + np.piecewise(t,[t<=pm,(pm<t)*(t<p2),t>=p2],[0.0,-1.0,0.0])   

    def qs_asymmetric(self,t,n):
        p1 = self.ps + n/self.ftrain
        pm = p1 + self.pw * 0.5
        pe = p1 + self.pw * 2.5
        return np.piecewise(t,[t<=p1,(p1<t)*(t<pm),t>=pm],[0.0,+1.0,0.0]) + 0.25 * np.piecewise(t,[t<=pm,(pm<t)*(t<pe),t>=pe],[0.0,-1.0,0.0])  
    
    def qs_shifted(self,t,n):
        p1  = self.ps + n/self.ftrain
        pe1 = p1 + self.pw * 0.5
        pe2 = p1 + 0.5/self.ftrain
        pe3 = pe2 + 0.5*self.pw
        return np.piecewise(t,[t<=p1,(p1<t)*(t<pe1),t>=pe1],[0.0,+1.0,0.0]) + np.piecewise(t,[t<=pe2,(pe2<t)*(t<pe3),t>=pe3],[0.0,-1.0,0.0])

    # Given the above waveform, evaluate the function
    def qs_evaluate(self, t, qsfun):
        n = 0
        pot_qs = np.zeros(len(t), dtype = float)
        while n < self.n_shf:
            pot_qs += qsfun(t, n=n)
            n   +=1      
        return pot_qs 
