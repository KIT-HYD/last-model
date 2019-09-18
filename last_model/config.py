from ._utils import ConfigBaseMixin
from last_model import ext


class Config(ConfigBaseMixin):
    identifier = 'model_configuration'

    def check_validity(self):
        """Check Configuration validity

        Method checks if all needed configuration is set 
        and valid. This includes the four extension class 
        arrays are present and each name defined in them can
        be mapped to an extension.
        """
        mapping = dict()
        # check each array
        for array_name in ('setup', 'pre_main', 'post_main', 'output'):
                
            # map each name to an Extension class
            for extension in self.get(array_name):
                Cls = self._map_name(extension)

                 # store in internal mapping
                if extension not in mapping.keys():
                    mapping[extension] = Cls

        return mapping

    def _map_name(self, classname):
        # check if loading extension names from Global namespace is supported
        try:
            load_global = self.get('can_load_from_global')
        except KeyError:
            load_global = False

        # find the classname in ext module
        if hasattr(ext, classname):
            Cls = getattr(ext, classname)
        elif load_global:
            if classname in dir():
                Cls = eval(classname)
            else:
                raise RuntimeError('The Extension %s is not in ext and not in global.' % classname)
        else:
            # TODO: define a specific error here
            raise RuntimeError('The Extension %s is not known' % classname)

        # now Cls is the class and can be returned
        return Cls

    def get_extension_classes(self, Last):
        # get the mapping to class names
        mapping = self.check_validity()
        armed_mapping = dict()

        # initialize each mapping
        for key, ExtensionCls in mapping.items():
            armed_mapping[key] = ExtensionCls(Last)

        # set setup classes
        Last._setup_classes = [armed_mapping[key] for key in self.get('setup')]

        # set main classes
        Last._pre_main_classes = [armed_mapping[key] for key in self.get('pre_main')]
        Last._post_main_classes = [armed_mapping[key] for key in self.get('post_main')]
        Last._output_classes = [armed_mapping[key] for key in self.get('output')]
        
        # dev
        # -------------------------------------
        #infiltration = ext.Infiltration(Last)
        #prefflow = ext.PreferentialFlow(Last)
        #Last._setup_classes = [infiltration]
        #Last._pre_main_classes = [prefflow, infiltration]
        # -------------------------------------
