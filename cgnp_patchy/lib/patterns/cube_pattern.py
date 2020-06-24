from __future__ import division

import mbuild as mb
import numpy as np


class CubePattern(mb.Pattern):
    """A nanoparticle coating pattern where points are removed from points on six axies.

    Parameters
    ----------
    chain_density : float
        Density of chain coating on the nanoparticle (chains / nm^2)
    radius : float
        Radius of the nanoparticle (nm)
    fractional_sa : float
        Fractional surface area of the nanoparticle to exclude coating (nm^2)
    
    Note
    ----------
    - If fractional surface area is too large for the sphere you are trying to pattern,
        this program will not find any points that satisfy the cubic pattern, which will
        lead to errors. We've tested up to 0.8, which failed.
    - The issue happens when cutoff is close to 1
    """
    def __init__(self, chain_density, radius, fractional_sa, **args):
        if fractional_sa >= 0.8:
            raise Exception("Coating pattern 'cubic' only works for fraction surface area values of 0.8 and below.")

        pattern = mb.SpherePattern(int(chain_density * 4.0 * np.pi * radius**2.0))
        pattern.scale(radius)
        total_sa = 4.0 * np.pi * radius**2.0
        patch_sa = total_sa * fractional_sa
        cutoff = patch_sa / (8 * np.pi * radius)
        points = np.array([xyz for xyz in pattern.points if xyz[2] < radius-cutoff
                           and xyz[2] > cutoff-radius and xyz[1] < radius-cutoff
                           and xyz[1] > cutoff-radius and xyz[0] < radius-cutoff
                           and xyz[0] > cutoff-radius])
        
        super(CubePattern, self).__init__(points=points, orientations=None)

if __name__ == "__main__":
    from save_pattern import save_pattern
    cube_pattern = CubePattern(4.0, 2.5, 0.79)
    save_pattern('test.xyz', cube_pattern)
