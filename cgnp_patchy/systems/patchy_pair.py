import mbuild as mb
import numpy as np

class PatchyPair(mb.Compound):
    def __init__(self, nano, sep=4):
        super(PatchyPair, self).__init__()

        nano.translate_to(np.zeros(3))

        nano2 = mb.clone(nano)
        nano2.translate([sep, 0, 0])
        self.add(nano,'nano')
        self.add(nano2,'nano2')

        boundingbox = self.boundingbox
        boundingbox.lengths = boundingbox.lengths * 3
        self.periodicity = boundingbox.lengths

if __name__ == "__main__":
    from cgnp_patchy.cgnp_patchy import cgnp_patchy
    from cgnp_patchy.lib.chains import CGAlkane
    chain_proto = CGAlkane(n=6, cap_front=False, cap_end=True)
    tnp = cgnp_patchy(radius=2.5, bead_diameter=0.6, chain=chain_proto, chain_density=2.5, coating_pattern='random', seed=57)
    patchy_pair = PatchyPair(tnp)
    patchy_pair.save('patchy_pair.mol2', overwrite=True)
