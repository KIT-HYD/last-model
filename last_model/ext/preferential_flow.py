import numpy as np

from last_model._utils import ExtensionBase

class PreferentialFlow(ExtensionBase):
    """Basic Preferential Flow Extension

    Assumes preferential flow domain with certain number 
    of macropores shaped like a straight circular cylinder. 
    Water particles are spherical shaped in a cubic storage.

    """
    identifier = 'Preferential Flow'

    def init_last(self):
        # TODO hier noch auseinanderhalten, was init und was setup ist
        pass

    def setup(self):
        # load preferential flow specific parameters
        params = self.params.get('preferential_flow')

        # maximum length of a macropore
        self.last.pfd_z = np.arange(0, params['depth'], params['depth_step'])

        # macropore grid elements
        self.last.pfd_dz = np.abs(np.diff(self.last.pfd_z))
        self.last.pfd_dim = self.last.pfd_z.size

        # macropore diameter
        self.last.pfd_D = params['diameter_mm'] / (10 * 100.)

        # amount of macropores within the domain
        self.last.n_mak = params['amount']

        # macropore radius
        self.last.pfd_r = self.last.pfd_D / 2. 

        # total particles within macropore
        self.last.pfd_N = params['total_particles']
        self.last.rho = 1000 # parameter?

        # initial soil moisture and solute concentration arrays
        # TODO: hier sollten wir array formate vereinheitlichen!
        self.last.pfd_theta = np.zeros((self.last.pfd_dim, 1))
        self.last.pfd_Cw_initial = np.zeros((self.last.pfd_dim, 1))

        # TODO: this is taken from the init_LAST 
        self.last.pfd_Cw = self.last.pfd_Cw_initial

        #-----------------------------------
        # Initialize macropore distribution
        #-----------------------------------

        # TODO: can this be put into a organized structure? dict?
        # TODO: I don't really get the parameters here...
        # depth distribution of macropores within pfd
        self.last.mak_big = self.last.pfd_z[self.last.pfd_dim - 1]
        
        # TODO: maybe implement this with default parameters?
        self.last.mak_mid = params['mak_mid']
        self.last.mak_sml = params['mak_sml']
        self.last.rate_big = params['rate_big']
        self.last.rate_mid = params['rate_mid']
        self.last.rate_sml = params['rate_sml']

        # -----------------
        # TODO: I think everything up to here was initialization, now follows setup
        # main difference: Other Extensions could alter initialized parameters...
        # TODO: some stuff needs to be rephrased.
        # TODO: I need a statement what of the following stuff needs to be attached
        # to LAST and what should be attached to PreferentialFlow

        # volume of grid element
        self.last.pfd_maxV_grid = np.pi * (self.last.pfd_r**2) * self.last.pfd_dz

        # total volume of a macropore
        self.last.V_gesamt = self.last.pfd_maxV_grid.sum()

        # water mass fitting in a grid element
        self.last.M_eventgrid = self.last.rho * self.last.pfd_maxV_grid

        # total water mass fitting in a macropore  
        self.last.m_gesamt = self.last.M_eventgrid.sum() 

        # mass of one particle
        self.last.pfd_m = self.last.m_gesamt / self.last.pfd_N

        # volume of a particle 
        self.last.V_particle = self.last.pfd_m / self.last.rho

        # diameter of a particle 
        self.last.pfd_D_particle = (self.last.V_particle / (np.pi/6))**(1/3) 

        # radius of a macropore
        self.last.r_particle = self.last.pfd_D_particle / 2. 

        # amount of particles fitting into a grid element
        self.last.pfd_n = self.last.M_eventgrid / self.last.pfd_m 

        # amount of particles fitting into a macropore
        self.last.n_gesamt = self.last.pfd_n.sum() 

        # circumference of a grid element
        self.last.U_eventgrid = 2 * np.pi * self.last.pfd_r 

        # amount of particles fitting with their diameter next to each other in a row 
        self.last.particle_contact = np.floor(self.last.U_eventgrid / self.last.pfd_D_particle)

        # amount of rows fitting over each other  
        self.last.anzahl_reihen = np.floor(self.last.pfd_dz[0] / self.last.pfd_D_particle) 

        # amount of particles having contact to lateral surface of grid element
        self.last.particle_contact_grid = self.last.particle_contact * self.last.anzahl_reihen 

        # total amount of particles having contact to lateral surface of a macropore
        self.last.particle_contact_total = self.last.particle_contact_grid * self.last.pfd_dz.size 

        # ensures that the amount of contact particles is not higher than the total amount of particles within the pfd
        if (self.last.particle_contact_total > self.last.pfd_N): 
            self.last.particle_contact_total = self.last.pfd_N
            self.last.particle_contact_grid = self.last.particle_contact_total / self.last.pfd_dz.size

        # flux density of a macropore
        # TODO: parameter? where does the 2884.2 come from?
        self.last.pfd_qmak = 2884.2 *(self.last.pfd_r**2) 
        
        # advective velocity of a particle 
        self.last.pfd_v = self.last.pfd_qmak 
        
        # hydraulic conductivity
        # TODO: why is this duplicated?
        self.last.pfd_k = self.last.pfd_v
        
        # no diffusion within macropore  
        # TODO: parameter?
        self.last.pfd_Diff = 0 
