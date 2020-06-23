from __future__ import division

import mbuild as mb
import numpy as np

def cartesian_to_spherical(origin=None, pos=None):
    """ Converts cartesian coordinates to spherical coordinates.

    Parameters
    ----------
    origin : np array
        Center of the sphere
    pos : np array
        Coordinates to convert
    """
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
 
    '''
    if (phi < 0):
        theta = theta + np.pi
        phi = abs(phi)
 
    if (num == 2):
        theta = (2*np.pi)-theta
    elif (num == 3):
        theta = (np.pi)+theta
    elif (num == 1):
        theta = (np.pi/2)+theta
    '''

    #theta = theta * (180/np.pi)
    #phi = phi * (180/np.pi)

    return np.array((r, theta, phi))

def spherical_to_cartesian(pos=None):
    """Converts spherical coordinates to cartesian coordinates.

    Parameters
    ----------
    pos : np array
        Coordinates to convert
    """
    pos = pos.reshape((3,))

    x = pos[0] * np.sin(pos[2]) * np.cos(pos[1])
    y = pos[0] * np.sin(pos[2]) * np.sin(pos[1])
    z = pos[0] * np.cos(pos[2])

    return np.array((x, y, z))

class TetrahedralPattern(mb.Pattern):
    """A nanoparticle coating pattern where points are removed from four poles.

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
        patch_cutoff = np.sqrt((patch_sa)/(4*np.pi)) 
        #cutoff = patch_sa / (8 * np.pi * radius)
        cutoff = patch_cutoff 

        #spherical_points = cartesian_to_spherical(pos=pattern.points)
        
        top_patch = np.array((radius, 0, 0))
        bottom_patch1 = np.array((radius, 0, (top_patch[2] + (109.5 * np.pi / 180))))
        bottom_patch2 = np.array((radius, (120 * np.pi / 180), bottom_patch1[2]))
        bottom_patch3 = np.array((radius, bottom_patch2[1] + (120 * np.pi / 180), bottom_patch1[2]))

        top_patch_cartesian = spherical_to_cartesian(pos=top_patch)
        bottom_patch1_cartesian = spherical_to_cartesian(pos=bottom_patch1)
        bottom_patch2_cartesian = spherical_to_cartesian(pos=bottom_patch2)
        bottom_patch3_cartesian = spherical_to_cartesian(pos=bottom_patch3)
 
        points = []

        for xyz in pattern.points:
            if xyz[0] > top_patch_cartesian[0]-patch_cutoff and xyz[0] < top_patch_cartesian[0]+patch_cutoff and xyz[1] > top_patch_cartesian[1]-patch_cutoff and xyz[1] < top_patch_cartesian[1]+patch_cutoff and xyz[2] > top_patch_cartesian[2]-patch_cutoff and xyz[2] < top_patch_cartesian[2]+patch_cutoff:
                points.append(xyz)

        for xyz in pattern.points:
            if xyz[0] > bottom_patch1_cartesian[0]-patch_cutoff and xyz[0] < bottom_patch1_cartesian[0]+patch_cutoff and xyz[1] > bottom_patch1_cartesian[1]-patch_cutoff and xyz[1] < bottom_patch1_cartesian[1]+patch_cutoff and xyz[2] > bottom_patch1_cartesian[2]-patch_cutoff and xyz[2] < bottom_patch1_cartesian[2]+patch_cutoff:
                points.append(xyz)

        for xyz in pattern.points:
            if xyz[0] > bottom_patch2_cartesian[0]-patch_cutoff and xyz[0] < bottom_patch2_cartesian[0]+patch_cutoff and xyz[1] > bottom_patch2_cartesian[1]-patch_cutoff and xyz[1] < bottom_patch2_cartesian[1]+patch_cutoff and xyz[2] > bottom_patch2_cartesian[2]-patch_cutoff and xyz[2] < bottom_patch2_cartesian[2]+patch_cutoff:
                points.append(xyz)

        for xyz in pattern.points:
            if xyz[0] > bottom_patch3_cartesian[0]-patch_cutoff and xyz[0] < bottom_patch3_cartesian[0]+patch_cutoff and xyz[1] > bottom_patch3_cartesian[1]-patch_cutoff and xyz[1] < bottom_patch3_cartesian[1]+patch_cutoff and xyz[2] > bottom_patch3_cartesian[2]-patch_cutoff and xyz[2] < bottom_patch3_cartesian[2]+patch_cutoff:
                points.append(xyz)


        '''for xyz in pattern.points:
            if xyz[0] < (cutoff-radius)-top_patch_cartesian[0] and xyz[1] < (cutoff-radius)-top_patch_cartesian[1] and xyz[2] < (cutoff-radius)-top_patch_cartesian[2]:
                points.append(xyz)
            else:
                continue


        for xyz in pattern.points:
            if xyz[0] > (cutoff-radius)-top_patch_cartesian[0]: # or xyz[0] > cutoff-top_patch_cartesian[0]:
                if xyz[1] > (cutoff-radius)-top_patch_cartesian[1]: # or xyz[1] > cutoff-top_patch_cartesian[1]:
                    if xyz[2] > (cutoff-radius)-top_patch_cartesian[2]: # or xyz[2] > cutoff-top_patch_cartesian[2]:
                        points.append(xyz)
                    else:
                        continue
                else:
                    continue
            else:
                continue


        for xyz in pattern.points:
            if xyz[0] < (cutoff)-bottom_patch1_cartesian[0] and xyz[1] < (cutoff)-bottom_patch1_cartesian[1] and xyz[2] < (cutoff)-bottom_patch1_cartesian[2]:
                points.append(xyz)
            else:
                continue
        
        for xyz in pattern.points:
            if xyz[0] < (cutoff-radius)-bottom_patch1_cartesian[0]:# or xyz[0] > cutoff-bottom_patch1_cartesian[0]:
                if xyz[1] < (cutoff-radius)-bottom_patch1_cartesian[1]: # or xyz[1] > cutoff-bottom_patch1_cartesian[1]:
                    if xyz[2] < (cutoff-radius)-bottom_patch1_cartesian[2]: # or xyz[2] > cutoff-bottom_patch1_cartesian[2]:
                        points.append(xyz)
                    else:
                        continue
                else:
                    continue
            else:
                continue

        for xyz in pattern.points:
            if xyz[0] > (cutoff-radius)-bottom_patch1_cartesian[0]:# or xyz[0] > cutoff-bottom_patch1_cartesian[0]:
                if xyz[1] > (cutoff-radius)-bottom_patch1_cartesian[1]: # or xyz[1] > cutoff-bottom_patch1_cartesian[1]:
                    if xyz[2] > (cutoff-radius)-bottom_patch1_cartesian[2]: # or xyz[2] > cutoff-bottom_patch1_cartesian[2]:
                        points.append(xyz)
                    else:
                        continue
                else:
                    continue
            else:
                continue

        for xyz in pattern.points:
            if xyz[0] < (cutoff)-bottom_patch2_cartesian[0] and xyz[1] < (cutoff)-bottom_patch2_cartesian[1] and xyz[2] < (cutoff)-bottom_patch2_cartesian[2]:
                points.append(xyz)
            else:
                continue

        for xyz in pattern.points:
            if xyz[0] < (cutoff-radius)-bottom_patch2_cartesian[0]: # or xyz[0] > cutoff-bottom_patch2_cartesian[0]:
                if xyz[1] < (cutoff-radius)-bottom_patch2_cartesian[1]: # or xyz[1] > cutoff-bottom_patch2_cartesian[1]:
                    if xyz[2] < (cutoff-radius)-bottom_patch2_cartesian[2]: # or xyz[2] > cutoff-bottom_patch2_cartesian[2]:
                        points.append(xyz)
                    else:
                        continue
                else:
                    continue
            else:
                continue

        for xyz in pattern.points:
            if xyz[0] > (cutoff-radius)-bottom_patch2_cartesian[0]:# or xyz[0] > cutoff-bottom_patch1_cartesian[0]:
                if xyz[1] > (cutoff-radius)-bottom_patch2_cartesian[1]: # or xyz[1] > cutoff-bottom_patch1_cartesian[1]:
                    if xyz[2] > (cutoff-radius)-bottom_patch2_cartesian[2]: # or xyz[2] > cutoff-bottom_patch1_cartesian[2]:
                        points.append(xyz)
                    else:
                        continue
                else:
                    continue
            else:
                continue

        for xyz in pattern.points:
            if xyz[0] > (cutoff-radius)-bottom_patch3_cartesian[0]: # or xyz[0] > cutoff-bottom_patch3_cartesian[0]:
                if xyz[1] > (cutoff-radius)-bottom_patch3_cartesian[1]: # or xyz[1] > cutoff-bottom_patch3_cartesian[1]:
                    if xyz[2] > (cutoff-radius)-bottom_patch3_cartesian[2]: # or xyz[2] > cutoff-bottom_patch3_cartesian[2]:
                        points.append(xyz)
                    else:
                        continue
                else:
                    continue
            else:
                continue

        for xyz in pattern.points:
            if xyz[0] < (cutoff-radius)-bottom_patch3_cartesian[0]: # or xyz[0] > cutoff-bottom_patch3_cartesian[0]:
                if xyz[1] < (cutoff-radius)-bottom_patch3_cartesian[1]: # or xyz[1] > cutoff-bottom_patch3_cartesian[1]:
                    if xyz[2] < (cutoff-radius)-bottom_patch3_cartesian[2]: # or xyz[2] > cutoff-bottom_patch3_cartesian[2]:
                        points.append(xyz)
                    else:
                        continue
                else:
                    continue
            else:
                continue
        
        for xyz in pattern.points:
            if xyz[0] < (cutoff)-bottom_patch3_cartesian[0] and xyz[1] < (cutoff)-bottom_patch3_cartesian[1] and xyz[2] < (cutoff)-bottom_patch3_cartesian[2]:
                points.append(xyz)
            else:
                continue'''
       
        # This pattern was written backwards, meaning that it returned the patches instead of everything but the patches. The code below is a temporary way to fix this until the math can be altered.
        reverse_pattern = []
        for point in pattern:
            if not np.all(np.isin(point, points)):
                reverse_pattern.append(point)

        super(TetrahedralPattern, self).__init__(points=reverse_pattern, orientations=None)


if __name__ == "__main__":
    from save_pattern import save_pattern
    #test = cartesian_to_spherical(pos=np.array((6, 2, -8)))
    tetrahedral_pattern = TetrahedralPattern(3.5, 3.0, 0.1)
    save_pattern('test_tetrahedral.xyz', tetrahedral_pattern)
