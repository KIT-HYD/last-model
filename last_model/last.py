from datetime import datetime as dt

from .config import Config
from .params import Params


class Last:
    def __init__(self, config=None, params=None):
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
        if config is None:
            self._config = Config()
        else:
            self._config = config

        # Set the Params module
        if params is None:
            self.params = Params()
        else:
            self.params = params
        
        # setup the Extension arrays
        self._setup_classes = []
        self._pre_main_classes = []
        self._post_main_classes = []
        self._output_classes = []

        # workflow
        self.workflow = []

        # fill the Extension arrays
        self._config.get_extension_classes(self)

        # run the init functions
        for instance in set([*self._setup_classes, *self._pre_main_classes, *self._post_main_classes, *self._output_classes]):
            instance._init_last()

        # PARAMETERS
        #-----------
        # break condition
        self.t_end = self.params.get('t_end')
        self.dtc = self.params.get('dtc')
        self.time = 0

        # profiling - as of now just runtime
        self.instance_creation = dt.now()
        self.main_start = None
        self.main_end = None
        self.did_run = False
    
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
            self.workflow.append('[%s] MAIN ALGORITHM' % str(dt.now()))
            self.main()

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
