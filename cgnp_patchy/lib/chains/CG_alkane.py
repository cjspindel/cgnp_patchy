import mbuild as mb
import numpy as np
from cgnp_patchy.lib.moieties import MMM, MME

class CG_alkane(mb.Compound):
    """ A coarse-grained, alkane chain, with three carbons to a CG bead, optionally terminated by CH3 particles """
    def __init__(self, n=6, cap_front=False, cap_end=True):
        """ Initalize a CG_alkane Compound.

        Parameters:
        ----------
        n : int
            Number of cg alkane particles in chain
        cap_front : boolean
            Add methyl group to beginning of chain ('down' port)
        cap_end : boolean
            Add methyl group to end of chain ('up' port)
        """

        super(CG_alkane, self).__init__()

        # Adjust length of Polymer for absense of methyl terminations.
        if not cap_front:
            n += 1
        if not cap_end:
            n += 1
        chain = mb.recipes.Polymer(MMM(), n=n-2, port_labels=('up', 'down'))
        self.add(chain, 'chain')
        if cap_front:
            self.add(MME(), 'methyl_front')
            mb.force_overlap(self['chain'], self['chain']['up'],
                             self['methyl_front']['up'])
        else:
            # Hoist port label to CG_alkane level.
            self.add(chain['up'], label='up', containment=False)

        if cap_end:
            self.add(MME(), 'methyl_end')
            mb.force_overlap(self['methyl_end'], self['methyl_end']['up'],
                             self['chain']['down'])
        else:
            # Hoist port label to CG_alkane level.
            self.add(chain['down'], label='down', containment=False)

