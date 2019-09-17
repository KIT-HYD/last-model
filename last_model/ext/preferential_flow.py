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

        # initial soil moisture and solute concentration arrays
        # TODO: hier sollten wir array formate vereinheitlichen!
        self.last.pfd_theta = np.zeros((self.last.pfd_dim, 1))
        self.last.pfd_Cw_initial = np.zeros((self.last.pfd_dim, 1))

        # TODO weiter in prefflow_domain.py ab Zeile 31

