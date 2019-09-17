from datetime import datetime as dt

import numpy as np

from .config import Config
from .params import Params
from .data_loader import DataLoader


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

        # TODO weiter in init_LAST.py Zeile 46

        #--------------------------------------
        # READ DATA
        #--------------------------------------
        self.read_data()

        # EXTENSION initialization
        # As a last step - run the init functions
        for instance in set([*self._setup_classes, *self._pre_main_classes, *self._post_main_classes, *self._output_classes]):
            instance._init_last()

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
