from last_model._utils import ExtensionBase


class Infiltration(ExtensionBase):
    def __init__(self, Last):
        super(Infiltration, self).__init__(Last)
        self.id = 'Infiltration'
