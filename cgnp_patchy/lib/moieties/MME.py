import mbuild as mb
import numpy as np

class MME(mb.Compound):
    """ Coarse-grained alkane bead containing a CH2-CH2-CH3 group """
    def __init__(self):
        super(MME, self).__init__()
        self.add(mb.Particle(name="_MME"))
        self.add(mb.Port(anchor=self[0], orientation=[0, 1, 0], separation=0.15), 'up')

