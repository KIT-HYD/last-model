import os

import numpy as np
import pandas as pd

from ._utils import ConfigBaseMixin
from .read_functions import (
    read_precip,
    read_precip_conc,
    read_init_theta,
    read_init_Cw
)


class DataLoader(ConfigBaseMixin):
    identifier = 'data'

    def get_data(self, key):
        """Load data

        Load the data as specified in the configuration file.
        The Module will automatically detect absolute and 
        relative paths to data files. Alternatively, the data can
        directly be given in the config file. 

        For the standard precipitation and isotope data needed for 
        the core LAST model, the input data can also be constructed 
        from given data.

        """
        # get the value
        value = self.get(key)

        if isinstance(value, str):
            abspath = os.path.abspath(os.path.join(os.path.dirname(self.path), value))

            # check if the string is a relative path
            if os.path.exists(abspath):
                return self._read_from_file(key, abspath)

            # check if the string is a absolute path
            elif os.path.exists(value):
                return self._read_from_file(key, value)

            else:
                raise NotImplementedError('Only absoult and relative paths are supported.')

        elif isinstance(value, list):
            # create array
            array = np.asarray(value)
            return self._read_from_array(key, array)

        elif isinstance(value, dict):
            raise NotImplementedError

        else:
            raise ValueError('Cannot parse the %s value' % key)

    def _read_from_file(self, key, path):
        # read data
        data = pd.read_csv(path, sep=' ', header=None).values
        
        # return using array function
        return self._read_from_array(key=key, array=data)


    def _read_from_array(self, key, array):
        if key.lower() == 'precipitation':
            return read_precip(array)
        elif key.lower() == 'precipitation_concentration':
            return read_precip_conc(array)
        elif key.lower() == 'soil_moisture':
            return read_init_theta(array, z=self.last.z, dim=self.last.dim)
        elif key.lower() == 'concentration':
            return read_init_Cw(array, z=self.last.z, dim=self.last.dim)
        else:
            return array
