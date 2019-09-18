import os
import json
import uuid
from datetime import datetime as dt


class ConfigBaseMixin:
    identifier = None

    def __init__(self, Last, path=None):
        if self.identifier is None:
            raise RuntimeError('Do not use ConfBaseMixin directly. Overwrite identifier, to mark file keyword to use.')
        
        # Last model reference
        self.last = Last

        # if no path is given, the current directory is assumed
        if path is None:
            path = os.getcwd()
        
        # append the parameter file
        if not path.endswith('.json'):
            path = os.path.join(path, 'config.json')
        
        if not os.path.exists(path):
            raise ValueError('The config file was not found at %s.' % os.path.dirname(path))
        else:
            # save path
            self.path = path
            self._config = dict()
        
        # load the file content
        self.load()

    def set(self, key, value):
        self._config[key] = value

    def get(self, key):
        return self._config[key]

    def _read_file(self):
        with open(self.path, 'r') as fs:
            return json.load(fs)

    def load(self):
        # load file content
        conf = self._read_file()
        
        params = conf.get(self.identifier)

        # set as parameters
        for key, value in params.items():
            self.set(key, value)
        
    def save(self):
 
        conf = self._read_file()

        # set the new parameters
        conf[self.identifier] = self._params

        # save _params to file
        with open(self.path, 'w') as fs:
            json.dump(fs, conf, indent=4)


class ExtensionBase:
    """Extension Base class

    This class can be used to implement custom algorithm parts.
    It accepts the LAST model instance as first argument and 
    makes it accessable as last attribute. Additionally, 
    the associated params object is made available as params. 
    Each Extension needs a init_last, setup and run method. 
    These methods will be called by LAST according to the
    configuration done by the user in the appropiate places.

    """
    identifier = None

    def __init__(self, Last):
        # bind the model
        self.last = Last
        self.params = Last.params  
        
        # set a unique identifier, if still None
        if self.identifier is None:
            self.identifier = uuid.uuid4()

    def _init_last(self):
        msg = '[%s] %s init_last function.'
        
        # run the user defined function
        try:
            self.init_last()
        except AttributeError:
            msg += ' Not defined.'
        
        # append to workflow
        self.last.workflow.append(msg % (str(dt.now()), self.identifier))

    def _setup(self):
        msg = '[%s] %s setup function.'

        # run the user defined function
        try:
            self.setup()
        except AttributeError as e:
            if "has no attribute 'setup'" in str(e):
                msg += ' Not defined.'
            else:
                raise e
        
        # append to workflow
        self.last.workflow.append(msg % (str(dt.now()), self.identifier))

    def _run(self):
        msg = '[%s] %s run function.'

        # run the user defined function
        try:
            self.run()
        except AttributeError:
            msg += ' Not definded.'
        
        # append to workflow
        self.last.workflow.append(msg % (str(dt.now()), self.identifier))