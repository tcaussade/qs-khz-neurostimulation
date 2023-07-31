import numpy as np


class Tissue:

    def __init__(self, tissue_type):
        self.eps_inf      = self.select_tissue(tissue_type)["einf"]
        self.sigma_ionic  = self.select_tissue(tissue_type)["s0"]
        self.eps          = self.select_tissue(tissue_type)["eps"]
        self.tau          = self.select_tissue(tissue_type)["tau"]
        self.alpha        = self.select_tissue(tissue_type)["alpha"]

        self.eps0 = 8.854187817e-12
        self.mu0  = 4*np.pi*1e-7

        pass

    #########################
    # Dispersion computations
    #########################

    def cole_cole_model(self, w):     
        # w must be angular frequency
        epsc = self.eps_inf \
                + np.sum(self.eps/(1 + (1j*w*self.tau)**(1-self.alpha) )) \
                + self.sigma_ionic/(1j*w*self.eps0)
        # conductivity is in [S/m]
        return np.array([self.eps0*np.real(epsc), -w*self.eps0*np.imag(epsc)])


    def select_tissue(self,tissue):

        if tissue == "brain_grey":
            # Brain grey-matter parameters (according to Gabriel dispersion model)
            eps_inf     = 4.0
            sigma_ionic = 0.02
            eps         = np.array([4.5e1,4.0e2,2.0e5,4.5e7], dtype = float)
            tau         = np.array([7.958e-12, 1.592e-8, 1.061e-4, 5.305e-3], dtype = float)
            alpha       = np.array([0.1,0.15,0.22,0.0])

        elif tissue == "brain_white":
            # Brain white-matter parameters (according to Gabriel dispersion model)
            eps_inf     = 4.0
            sigma_ionic = 0.02
            eps         = np.array([32,100, 4.0e4,3.5e7], dtype = float)
            tau         = np.array([7.96e-12, 7.96e-9, 53.05e-6, 7.958e-3], dtype = float)
            alpha       = np.array([0.1,0.1,0.3,0.02], dtype  = float)
        
        elif tissue == "liver":
            # Liver parameters (according to Gabriel dispersion model)
            eps_inf     = 4.0
            sigma_ionic = 0.02
            eps         = np.array([39,6000, 5.0e4,3.0e7], dtype = float)
            tau         = np.array([8.84e-12, 530.52e-9, 22.74e-6, 15.915e-3], dtype = float)
            alpha       = np.array([0.1,0.2,0.2,0.05], dtype  = float)
        
        elif tissue == "bone":
            # cortical bONE parameters (according to Gabriel dispersion model)
            eps_inf     = 2.5
            sigma_ionic = 0.02
            eps         = np.array([10,180, 5.0e3,1.0e5], dtype = float)
            tau         = np.array([13.26e-12, 79.58e-9, 159.15e-6, 15.915e-3], dtype = float)
            alpha       = np.array([0.2,0.2,0.2,0.0], dtype  = float)

        elif tissue == "skin":
            # skin wet parameters
            eps_inf     = 4.0
            sigma_ionic = 0.0004
            eps         = np.array([39,280, 3e4,3e4], dtype = float)
            tau         = np.array([7.96e-12, 79.58e-9, 1.59e-6, 1.592e-3], dtype = float)
            alpha       = np.array([0.1,0.0,0.16,0.2], dtype  = float)

        elif tissue == "fat":
            # not infiltrated
            eps_inf     = 2.5
            sigma_ionic = 0.01
            eps         = np.array([3, 15, 3.3e4,1e7], dtype = float)
            tau         = np.array([7.96e-12, 15.92e-9, 159.15e-6, 7.958e-3], dtype = float)
            alpha       = np.array([0.2,0.1,0.05,0.01], dtype  = float)

        elif tissue == "muscle":
            # not infiltrated
            eps_inf     = 4.0
            sigma_ionic = 0.02
            eps         = np.array([50, 7000, 1.2e6,2.5e7], dtype = float)
            tau         = np.array([7.23e-12, 353.68e-9, 318.31e-6, 2.274e-3], dtype = float)
            alpha       = np.array([0.1,0.1,0.1,0.00], dtype  = float)

        else:
            print("NO TISSUE HAS BEEN SELECTED")
            print("\ncheck your spelling...")

        return {"einf": eps_inf,
                "s0":   sigma_ionic,
                "eps":  eps,
                "tau":  tau,
                "alpha":alpha}
