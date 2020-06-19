from __future__ import division

import mbuild as mb
import numpy as np


def cartesian_to_spherical(origin=None, pos=None):
    if origin==None:
        origin = np.array((0, 0, 0))
    else:
        origin = np.ndarray(origin)

    #pos = np.ndarray(pos, dtype=float)

    pos = pos.reshape((3,))    
    origin = origin.reshape((3,))

    r = np.linalg.norm(pos-origin)

    phi = np.arccos(pos[2]/r)
   
    theta = np.arctan(pos[1]/pos[0])

    num = sum(1 for number in pos if number < 0)
   
    if (pos[0] < 0 and pos[1] < 0):
        theta = theta - (np.pi) + (2*np.pi)
    elif (pos[0] < 0 and pos[1] > 0):
        theta = theta + (np.pi)
    elif (pos[1] < 0 and pos[0] > 0):
        theta = (2*np.pi + theta) 
    else:
        pass
 
    #theta = theta * (180/np.pi)
    #phi = phi * (180/np.pi)

    return np.array((r, theta, phi))

def spherical_to_cartesian(pos=None):

    pos = pos.reshape((3,))

    x = pos[0] * np.sin(pos[2]) * np.cos(pos[1])
    y = pos[0] * np.sin(pos[2]) * np.sin(pos[1])
    z = pos[0] * np.cos(pos[2])

    return np.array((x, y, z))

class RingPattern(mb.Pattern):
    """A pattern where points are removed from three poles. Tetrahedral pattern without top patch.

    Parameters
    ----------
    chain_density : float
        Density of chain coating on the nanoparticle (chains / nm^2)
    radius : float
        Radius of the nanoparticle (nm)
    fractional_sa : float
        Fractional surface area of the nanoparticle to exclue coating (nm^2)
    """
    def __init__(self, chain_density, radius, fractional_sa, **args):
        
        pattern = mb.SpherePattern(int(chain_density * 4.0 * np.pi * radius**2.0))
        pattern.scale(radius)
        total_sa = 4.0 * np.pi * radius**2.0
        patch_sa = total_sa * fractional_sa
        patch_cutoff = np.sqrt((patch_sa)/(4*np.pi)) 
        #cutoff = patch_sa / (8 * np.pi * radius)
        cutoff = patch_cutoff 
        print(total_sa)
        print(cutoff)
        print(patch_sa)

        #spherical_points = cartesian_to_spherical(pos=pattern.points)
        
        bottom_patch1 = np.array((radius, 0, (0 + (109.5 * np.pi / 180))))
        bottom_patch2 = np.array((radius, (120 * np.pi / 180), bottom_patch1[2]))
        bottom_patch3 = np.array((radius, bottom_patch2[1] + (120 * np.pi / 180), bottom_patch1[2]))

        bottom_patch1_cartesian = spherical_to_cartesian(pos=bottom_patch1)
        bottom_patch2_cartesian = spherical_to_cartesian(pos=bottom_patch2)
        bottom_patch3_cartesian = spherical_to_cartesian(pos=bottom_patch3)
       
        print(bottom_patch1_cartesian)
        print(bottom_patch2_cartesian)
        print(bottom_patch2_cartesian)
 
        points = []

        '''
        for xyz in pattern.points:
            if xyz[0] > top_patch_cartesian[0]-patch_cutoff and xyz[0] < top_patch_cartesian[0]+patch_cutoff and xyz[1] > top_patch_cartesian[1]-patch_cutoff and xyz[1] < top_patch_cartesian[1]+patch_cutoff and xyz[2] > top_patch_cartesian[2]-patch_cutoff and xyz[2] < top_patch_cartesian[2]+patch_cutoff:
                points.append(xyz)
        '''
        for xyz in pattern.points:
            if xyz[0] > bottom_patch1_cartesian[0]-patch_cutoff and xyz[0] < bottom_patch1_cartesian[0]+patch_cutoff and xyz[1] > bottom_patch1_cartesian[1]-patch_cutoff and xyz[1] < bottom_patch1_cartesian[1]+patch_cutoff and xyz[2] > bottom_patch1_cartesian[2]-patch_cutoff and xyz[2] < bottom_patch1_cartesian[2]+patch_cutoff:
                points.append(xyz)

        for xyz in pattern.points:
            if xyz[0] > bottom_patch2_cartesian[0]-patch_cutoff and xyz[0] < bottom_patch2_cartesian[0]+patch_cutoff and xyz[1] > bottom_patch2_cartesian[1]-patch_cutoff and xyz[1] < bottom_patch2_cartesian[1]+patch_cutoff and xyz[2] > bottom_patch2_cartesian[2]-patch_cutoff and xyz[2] < bottom_patch2_cartesian[2]+patch_cutoff:
                points.append(xyz)

        for xyz in pattern.points:
            if xyz[0] > bottom_patch3_cartesian[0]-patch_cutoff and xyz[0] < bottom_patch3_cartesian[0]+patch_cutoff and xyz[1] > bottom_patch3_cartesian[1]-patch_cutoff and xyz[1] < bottom_patch3_cartesian[1]+patch_cutoff and xyz[2] > bottom_patch3_cartesian[2]-patch_cutoff and xyz[2] < bottom_patch3_cartesian[2]+patch_cutoff:
                points.append(xyz)

        super(RingPattern, self).__init__(points=points, orientations=None)


if __name__ == "__main__":
    from save_pattern import save_pattern
    #test = cartesian_to_spherical(pos=np.array((6, 2, -8)))
    ring_pattern = RingPattern(3.5, 3.0, 0.1)
    save_pattern('test_ring.xyz', ring_pattern)


