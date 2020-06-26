import mbuild as mb
import numpy as np
from cgnp_patchy.lib.patterns import *
from cgnp_patchy.lib.utils.save_pattern import save_pattern

def count_patch_points(pattern, radius, chain_density):
    ''' Counts and returns the amount of points that are being removed in a certain nanoparticle coating pattern.

    Parameters
    ----------
    pattern : mb.Pattern
        Nanoparticle coating pattern to count which points were removed
    radius : float
        Radius of the nanoparticle
    chain_density : float
        Density of chains on the nanoparticle surface
    '''
    isotropic_pattern = mb.SpherePattern(int(chain_density*4.0*np.pi*radius**2))
    isotropic_pattern.scale(radius)
     
    patch = []
    for point in isotropic_pattern.points:
        if not np.all(np.isin(point, pattern.points)):
            patch.append(point)

    return len(patch)
