from last_model._utils import ExtensionBase

class PreferentialFlow(ExtensionBase):
    def __init__(self, Last):
        super(PreferentialFlow, self).__init__(Last)
        self.id = 'Preferential Flow'
