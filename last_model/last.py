from datetime import datetime as dt

import numpy as np
import pandas as pd

from .config import Config
from .params import Params
from .data_loader import DataLoader
from .init_functions import (
    psi_theta, 
    k_psi,
    init_particles,
    init_particles_pos
)


class Last:
    def __init__(self, path=None, **kwargs):
        """LAST model 

        Main class implementing LAST model. 
        You need two configuration objects, Config and Params
        Config will be used to set up the model logic, 
        Params define the model parameters to be used.

        Config
        ------

        Params
        ------

        Lifecycle
        ---------
        A call to LAST's run method will be processed in two main 
        phases. The Classes to be initialized and run are defined
        in the Config module.
        
        1. setup
        2. main

        The setup phase is run once to initialize border conditions,
        prepopulate data, warmup or reserve resources.
        The main phase will simulate until the break condition is met.

        """
        # Set the Config module
        self._config = Config(Last=self, path=path)

        # Set the Params module
        self.params = Params(Last=self, path=path)

        # Set the DataLoader module
        self.data = DataLoader(Last=self, path=path)
        
        # setup the Extension arrays
        self._setup_classes = []
        self._pre_main_classes = []
        self._post_main_classes = []
        self._output_classes = []

        # workflow
        self.workflow = []

        # fill the Extension arrays
        self._config.get_extension_classes(self)

        #--------------------------------------
        # PARAMETERS
        #--------------------------------------
        # break condition
        self.t_end = self.params.get('t_end')
        self.dtc = self.params.get('dtc')
        self.time = 0

        # profiling - as of now just runtime
        self.instance_creation = dt.now()
        self.main_start = None
        self.main_end = None
        self.did_run = False

        # soil matrix parameters
        self.z = np.arange(0, self.params.get('grid_depth'), self.params.get('grid_step'))
        self.dz = np.abs(np.diff(self.z))
        self.dim = self.z.size

        # soil parameters
        # TODO: die Namen sollten angepasst werden
        soil = self.params.get('soil')
        self.ks = soil.get('ks')
        self.ths = soil.get('ths')
        self.thr = soil.get('thr')
        self.alph = soil.get('alph')
        self.n_vg = soil.get('n_vg')
        self.stor = soil.get('stor')
        self.l_lg = soil.get('l_lg')
        
        # particle settings
        self.mob_fak = self.params.get('mobile_particle_fraction')
        self.N = self.params.get('total_particles')
        self.nclass = self.params.get('classes')

        # tmix distributiion
        self.tmix = self.params.get('tmix')

        # event particles
        # TODO: I think we can pre-allocate these arrays to make them faster

        # position of event particles entering the soil matrix
        self.position_z_event = np.array([])
        
        # age of event particles particles entering the soil matrix
        self.age_event = np.array([np.nan])

        # concentration of event particles entering the soil matrix
        self.Cw_event = np.array([np.nan])
        self.t_mix_event = np.array([np.nan])

        # TODO: Move the pfd_ into the PreferentialFlow extension?
        # position of particles within the pfd
        self.pfd_position_z = np.array([np.nan])

        # age of particles within the pfd
        self.pfd_age = np.array([np.nan])

        # concentration of particles entering pfd
        self.pfd_Cw_event = np.array([np.nan])

        # water surface storage
        self.m_surface = 0

        # solute surface storage
        self.m_slt_surface = 0
        self.m_input = 0    

        # TODO: can this be removed?
        # counter for amount of incoming precipitation water
        self.mass_totalwaterinput = 0

        # counter for amount of incoming solute mass
        self.mass_totalsoluteinput = 0

        # counter for water mass entering the soil matrix and the pfd
        # TODO: this should REALLY be renamed
        self.bla = 0 
        self.bla1 = 0 

        # counter for soulte mass entering the soil matrix and pfd
        self.solute_test = 0 
        self.solute_test2 = 0 

        self.theta_event = 0

        # TODO: this is calculated, maybe put it into an setup function
        # or into an extension
        # TODO: maybe only store and pass the soil dictionary?
        self.psi_init, self.c_init = psi_theta(alph=self.alph, dim=self.dim, n_vg=self.n_vg, stor=self.stor, theta=self.theta_init, thr=self.thr, ths=self.ths)

        # calculates initial hydraulic conductivities using psi
        self.k_init = k_psi(alph=self.alph, dim=self.dim, ks=self.ks, l_vg=self.l_vg, n_vg=self.n_vg, psi=self.psi_init)
        
        # create lookup table for D and k values
        self.prob = 1.0  # TODO: what is this? Parameter?
        self.D_table = np.zeros((1, self.nclass))
        self.K_table = np.zeros((1, self.nclass))

        # soil moisture bins
        self.theta_table = np.arange(self.thr, self.ths, ((self.ths - self.ths) / self.nclass)).transpose()

        # calculate initial psi from initial soil moisture
        for i in range(self.nclass):
            self._map_psi_to_theta(ibin=i)
        
        # now, the working parameters for each iteration can be initialized
        # TODO: do we need the *_init values for anything else than this step?
        # if no, we can move that part into a setup routine and remove them
        self.k = self.k_init
        self.c = self.c_init
        self.psi = self.psi_init
        self.Cw = self.Cw_init

        #-----------------------------------
        # PARTICLE INITIALIZATION
        #-----------------------------------
        # initialises particle mass (m) distribution (n), positions (position_z) and concentrations (c_particle)
        self.n, self.m = init_particles(theta=self.theta_init, dz=self.dz, dim=self.dim, N=self.N)

        # initialises particle positions and initial particle concentrations, particle ages, retarded and degraded particle solute concentration
        self.position_z, self.c_particle = init_particles_pos(z=self.z, dim=self.dim, n=self.n, Cw_init=self.Cw_init)
        self.age_particle = np.zeros((self.position_z.size, 1))

        # TODO: ui ui. Das machen wir anders
        self.position_z = pd.concat([pd.DataFrame(self.position_z), pd.DataFrame(self.age_particle), pd.DataFrame(c_particle)], axis=1)
        self.position_z = self.position_z.values

        # initial input conditions at soil surface
        self.Cw_eventparticles = self.prec_conc[0]

        # flux density of precipitation water input
        self.qb_u = -0.5 * (self.prec_int[0] + self.prec_int[1])

        # TODO plot settings are skipped for now...


        #--------------------------------------
        # READ DATA
        #--------------------------------------
        self.read_data()

        # time settings
        # TODO: here, I break with the original code, as the data 
        # is not loaded yet. I would not depend the break condition
        # on the input data...
        # on the long run, we should get rid of these settings and 
        # move them to parameters
        self.time = self.prec_time[0]
        self.t_end = self.t_end + self.time
        self.t_boundary = self.prec_time # TODO why exactly do we need it twice? ref. init_LAST.py lines 111 to 114
        self.i_time = 1




        # EXTENSION initialization
        # As a last step - run the init functions
        for instance in set([*self._setup_classes, *self._pre_main_classes, *self._post_main_classes, *self._output_classes]):
            instance._init_last()

    def _map_psi_to_theta(self, ibin):
        """
        TODO: for me it's not 100% clear what this function does
        """
        theta_actual = np.ones((self.dim, 1)) * self.theta_table[ibin]

        # calculate the psi of current bin?
        psi_h, c_h = psi_theta(alph=self.alph, dim=self.dim, n_vg=self.n_vg, stor=self.stor, theta=theta_actual, thr=self.thr, ths=self.ths)
        k_help = k_psi(alph=self.alph, dim=self.dim, ks=self.ks, l_vg=self.l_vg, n_vg=self.n_vg, psi=psi_h)

        # TODO: why is this always index 0 here? can we simplify this?
        if c_h[0] > 0:
            self.D_table[0][ibin] = k_help[0] / c_h[0]
        else:
            self.D_table[0][ibin] = 0
        self.K_table[0][i] = k_help[0]


    def read_data(self):
        """
        Load the minimal data needed for LAST to run 
        correctly
        """
        # precipitation
        self.prec_int, self.prec_time = self.data.get_data('precipitation')

        # concentration in precipitation
        self.prec_conc = self.data.get_data('precipitation_concentration')

        # initial soil moisture state
        self.theta_init = self.data.get_data('soil_moisture')

        # initial concentration in profile
        self.Cw_init = self.data.get_data('concentration')

        # final observed tracer concentration profile 
        self.concfin = self.data.get_data('concentration')
    
    def setup(self):
        """
        """
        for instance in self._setup_classes:
            instance._setup()
       
    def run(self):
        """
        """
        # setup the model
        self.setup()

        # run 
        self.main_start = dt.now()

        while self.time < self.t_end:
            self.time += self.dtc
            
            # run the pre main classes
            for instance in self._pre_main_classes:
                instance._run()
            
            # MAIN ALGORITHM
            self.main()
            self.workflow.append('[%s] MAIN ALGORITHM' % str(dt.now()))

            # run the post main classes
            for instance in self._post_main_classes:
                instance._run()

        # The model is finished.
        self.main_end = dt.now()
        self.did_run = True

        # run the output classes
        for instance in self._output_classes:
            instance._run()

    def main(self):
        pass
