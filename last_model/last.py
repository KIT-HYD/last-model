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
from .model_functions import (
    displ_mtx_pre,
    displ_mtx_event,
    displ_pfd,
    mixing_mtx,
    mixing_pfd_mtx
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

        # controls the displacement routine
        self.pc = None 

        # profiling - as of now just runtime
        self.instance_creation = dt.now()
        self.main_start = None
        self.main_end = None
        self.did_run = False

        # soil matrix parameters
        self.z = np.arange(0, self.params.get('grid_depth'), self.params.get('grid_step'))
        self.dz = np.abs(np.diff(self.z))
        self.dim = self.z.size

        # grid should be defined here, so load data
        # TODO: maybe not calculate the grid, but make it input data?
        #--------------------------------------
        # READ DATA
        #--------------------------------------
        self.read_data()

        # soil parameters
        # TODO: die Namen sollten angepasst werden
        soil = self.params.get('soil')
        self.ks = soil.get('ks')
        self.ths = soil.get('ths')
        self.thr = soil.get('thr')
        self.alph = soil.get('alph')
        self.n_vg = soil.get('n_vg')
        self.stor = soil.get('stor')
        self.l_vg = soil.get('l_vg')
        
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
        self.theta_table = np.arange(self.thr, self.ths, ((self.ths - self.thr) / self.nclass)).transpose()

        # calculate initial psi from initial soil moisture
        for i in range(self.nclass):
            self._map_psi_to_theta(ibin=i)
        
        # now, the working parameters for each iteration can be initialized
        # TODO: do we need the *_init values for anything else than this step?
        # if no, we can move that part into a setup routine and remove them
        self.k = self.k_init
        self.c = self.c_init
        self.psi = self.psi_init
        self.theta = self.theta_init
        self.Cw = self.Cw_init

        #-----------------------------------
        # PARTICLE INITIALIZATION
        #-----------------------------------
        # initialises particle mass (m) distribution (n), positions (position_z) and concentrations (c_particle)
        self.n, self.m = init_particles(theta=self.theta_init, dz=self.dz, dim=self.dim, N=self.N)

        # initialises particle positions and initial particle concentrations, particle ages, retarded and degraded particle solute concentration
        self.position_z, self.c_particle = init_particles_pos(z=self.z, dim=self.dim, n=self.n, Cw_init=self.Cw_init)
        
        # TODO: Why?
        self.position_z = np.round(self.position_z, 3)
        self.age_particle = np.zeros((self.position_z.size, 1))

        # TODO: ui ui. Das machen wir anders
        self.position_z = pd.concat([pd.DataFrame(self.position_z), pd.DataFrame(self.age_particle), pd.DataFrame(self.c_particle)], axis=1)
        self.position_z = self.position_z.values

        # initial input conditions at soil surface
        self.Cw_eventparticles = self.prec_conc[0]

        # flux density of precipitation water input
        self.qb_u = -0.5 * (self.prec_int[0] + self.prec_int[1])

        # TODO plot settings are skipped for now...

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
        self.K_table[0][ibin] = k_help[0]

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
        """Main Routine

        This function is run on each iteration until the break 
        condition is met. The Pre-Main classes are executed before
        this method is called and the Post-Main classes are 
        executed afterwards. In case particle dispalcement itself
        shall be intercepted, the follwoing lifecycle functions are
        available.

        The main difference between a pre-main class extension 
        and the pre-dispalcement hook is that the hook is 
        also run twice on each displacement and the pc parameter 
        is already set.

        Lifecycle:
        ----------

        1. Pre Particle displacement        [main_pre_displacement]
        2. Particle displacement            [main_displacement]
        3. Post particle displacement       [main_post_displacement]
        4. Mixing of matrix and macropores  [main_mixing]
        5. finalizing                       [main_finalizing]

        """
        self.pc = 0
        # run the displacement twice on half depth
        # TODO: is there a more elegant way?
        for pc in range(0,2):
            self.pc = pc

            # TODO: here, we could add profiling options in verbose mode
            # pre-displacement
            self.main_pre_displacement()

            # displacement routine
            self.main_displacement()

            # post displacement, with pc still in state
            self.main_post_displacement()

        # reset pc again
        self.pc = None

        # particle mixing
        self.main_mixing()

        # finish this iteration
        self.main_finalizing()

    def main_pre_displacement(self):
        """ 
        setup advective velocity and diffusivity from parameters
        then, set to 0 in each grid cell where the soil moisture 
        is close to thr (residual soil moisture)

        """
        # Step 1: initialisation of parameters v and D
        
        # advective velocity in each grid element of the soil matrix
        self.v = 1 * self.k
        # diffusivity in each grid element of the soil matrix
        self.D = 1 * self.k / self.c
        
        # finds all grid elements with a soil moisture near to thr (almost dry soil)
        # sets the velocity and diffusivity in this grid elements to 0-> no flux!
        ipres = self.theta < 1.1 * self.thr
        if ipres.any() == True:
            self.v[ipres==True] = 0
            self.D[ipres==True] = 0

    def main_displacement(self):
        """
        Displace the particles 
        """
        # displacement of pre-event particles in soil marix
        self.Cw, self.mtx_avg_age, self.position_z, self.theta = displ_mtx_pre(D=self.D, dim=self.dim, dtc=(0.5*self.dtc), D_table=self.D_table, dz=self.dz, K_table=self.K_table, m=self.m, mob_fak=self.mob_fak, position_z=self.position_z, ths=self.ths, v=self.v, z=self.z)

        # displacement of event particles in soil matrix
        if self.position_z_event.size > 0:
            self.position_z_event, self.theta_event = displ_mtx_event(dim=self.dim, dtc=(0.5*self.dtc), D_table=self.D_table, K_table=self.K_table, m=self.m, dz=self.dz, position_z_event=self.position_z_event, prob=self.prob, ths=ths, z=z)
        
        # displacement of particles in pfd and drainage; only one times in the predictor step because it is just dependent on constant advective veloscity
        if self.pc == 0 and np.isnan(self.pfd_position_z).all() == False:
            self.pfd_particles, self.pfd_Cw, self.pfd_theta, self.pfd_position_z = displ_pfd(n_mak=self.n_mak, pfd_dim=self.pfd_dim, pfd_dz=self.pfd_dz, pfd_m=self.pfd_m, pfd_n=self.pfd_n, pfd_position_z=self.pfd_position_z, pfd_r=self.pfd_r, pfd_z=self.pfd_z)       

    def main_post_displacement(self):
        """
        The psi, c and k parameters are recalculated after displacement
        """
        # update of parameters psi, c and k for the corrector step
        self.psi, self.c = psi_theta(alph=self.alph, dim=self.dim, n_vg=self.n_vg, stor=self.stor, theta=self.theta, thr=self.thr, ths=self.ths)
        self.k = k_psi(alph=self.alph, dim=self.dim, ks=self.ks, l_vg=self.l_vg, n_vg=self.n_vg, psi=self.psi)
 
    def main_mixing(self):
        """
        Mixing of event particle with pre-event particles 
        within the first grid elements of soil matrix
        """
        # TODO: when is that not true?
        if self.position_z_event.size > 0:
            self.age_event, self.Cw_event, self.position_z, self.position_z_event = mixing_mtx(age_event=self.age_event, Cw_event=self.Cw_event, position_z=self.position_z, position_z_event=self.position_z_event, time=self.time)

        val_pos = -np.sort(-self.position_z[:,0])
        idx = np.argsort(-self.position_z[:,0])    
        
        # TODO: I don't get this step. this should be made easier
        self.position_z = np.concatenate((np.reshape(val_pos, (val_pos[:,].size,1)),np.reshape(self.position_z[idx,1],(self.position_z[idx,1].size,1)),np.reshape(self.position_z[idx,2],(self.position_z[idx,2].size,1))), axis=1)               
    
        # Mixing between pfd and soil matrix
        if np.isnan(self.pfd_position_z).all() == False:
            # TODO: this has to be changed. parameters should be passed differently
            self.pfd_position_z, self.position_z = mixing_pfd_mtx(dtc=self.dtc, k=self.k, ks=self.ks, m=self.m, mak_sml=self.mak_sml, mak_mid=self.mak_mid, n_mak=self.n_mak, particle_contact_grid=self.particle_contact_grid, pfd_dim=self.pfd_dim, pfd_dz=self.pfd_dz, pfd_m=self.pfd_m, pfd_n=self.pfd_n, pfd_particles=self.pfd_particles, pfd_position_z=self.pfd_position_z, pfd_r=self.pfd_r, pfd_theta=self.pfd_theta, pfd_z=self.pfd_z, position_z=self.position_z, psi=self.psi, rate_big=self.rate_big, rate_mid=self.rate_mid, rate_sml=self.rate_sml, rho=self.rho, theta=self.theta, ths=self.ths, z=self.z)
 
    def main_finalizing(self):
        """
        Update input conditions at the top of the soil
        """
        # finds actual position in boundary conditions time series
        ipos = np.where(self.t_boundary <= self.time)[0][-1] 
        
        # updates flux density of precipitation water input
        self.qb_u = -0.5 * (self.prec_int[ipos,] + self.prec_int[ipos + 1,])
        
        # updates the concentration of precipitation water input 
        self.Cw_eventparticles = self.prec_conc[ipos,] 

