import mbuild as mb
import numpy as np
from cgnp_patchy.lib.moieties import MMM, MME

class CGAlkane(mb.Compound):
    def __init__(self, n=6, cap_front=False, cap_end=True):
        """ A coarse-grained, alkane chain with three carbons to a CG bead. Optionally terminated by CH3 particles.

        Parameters
        ----------
        n : int
            Number of cg alkane beads in chain. Beads are coarse-grained 3:1.
        cap_front : boolean
            Add methyl group to beginning of chain ('down' port). May currently cause problems attaching chain to nanoparticle.
        cap_end : boolean
            Add methyl group to end of chain ('up' port).
        """

        super(CGAlkane, self).__init__()

        # Adjust length of polymer for absense of methyl terminations.
        if not cap_front:
            n += 1
        if not cap_end:
            n += 1
        chain = mb.recipes.Polymer(MMM(), n=n-2, port_labels=('up', 'down'))
        self.add(chain, 'chain')
        
        # Add cap to front of chain. 
        if cap_front:
            self.add(MME(), 'methyl_front')
            mb.force_overlap(self['chain'], self['chain']['up'],
                             self['methyl_front']['up'])
        else:
            # Hoist port label to CGAlkane level.
            self.add(chain['up'], label='up', containment=False)

        # Add cap to end of chain.
        if cap_end:
            self.add(MME(), 'methyl_end')
            mb.force_overlap(self['methyl_end'], self['methyl_end']['up'],
                             self['chain']['down'])
        else:
            # Hoist port label to CGAlkane level.
            self.add(chain['down'], label='down', containment=False)

if __name__ == "__main__":
    chain = CGAlkane()
    chain.save('chain.mol2', overwrite=True)
