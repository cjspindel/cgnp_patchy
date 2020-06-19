from __future__ import division

import mbuild as mb
import numpy as np


class BipolarPattern(mb.Pattern):
    """A pattern where points are removed from two poles

    Parameters
    ----------
    chain_density : float
        Density of chain coating on the nanoparticle (chains / nm^2)
    radius : float
        Radius of the nanoparticle (nm)
    fractional_sa : float
        Fractional surface area of the nanoparticle to exclude coating (nm^2)
    """
    def __init__(self, chain_density, radius, fractional_sa, **args):
        pattern = mb.SpherePattern(int(chain_density * 4.0 * np.pi * radius**2.0))
        pattern.scale(radius)
        total_sa = 4.0 * np.pi * radius**2.0
        patch_sa = total_sa * fractional_sa
        cutoff = patch_sa / (4 * np.pi * radius)
        points = np.array([xyz for xyz in pattern.points if xyz[2] < radius-cutoff
                           and xyz[2] > cutoff-radius])

        super(BipolarPattern, self).__init__(points=points, orientations=None)

if __name__ == "__main__":
    from save_pattern import save_pattern
    bipolar_pattern = BipolarPattern(chain_density=4.0, radius=5.0, fractional_sa=0.5)
    save_pattern('test.xyz', bipolar_pattern)
