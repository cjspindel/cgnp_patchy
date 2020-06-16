from __future__ import division

import math

import mbuild as mb
import numpy as np
from scipy.spatial import distance

def _fast_sphere_pattern(n, radius):
    """Faster version of mBuild's SpherePattern. """
    phi = (1 + np.sqrt(5)) / 2
    long_incr = 2*np.pi / phi
    dz = 2.0 / float(n)
    bands = np.arange(n)
    z = bands * dz - 1.0 + (dz/2.0)
    r = np.sqrt(1.0 - z*z)
    az = bands * long_incr
    x = r * np.cos(az)
    y = r * np.sin(az)
    points = np.column_stack((x, y, z)) * np.asarray([radius])

    return points

class Nanoparticle(mb.Compound):
    """Coarse-gained nanoparticle class."""
    def __init__(self, r=5.0, sigma=0.8):
        super(Nanoparticle, self).__init__()

        r_CG = sigma / 2
        r_silica = 0.40323 / 2

        r = r - r_CG + r_silica

        # N_approx = a(R/sigma)^2 + b(R/sigma) + c
        # Not an exact number but very close
        a = 9.4379
        b = 0.6826
        c = -1.3333
        N_approx = a * ((r/sigma)**2) + b * (r/sigma) + c

        # Binary search algorithm to find maximum number of CG beads without overlaps
        min_points = max(N_approx - 500, 1)
        max_points = N_approx + 500
        opt_points = 0
        while opt_points == 0:
            mid = math.ceil((max_points + min_points) / 2)
            points = _fast_sphere_pattern(mid, r)
            points_high = _fast_sphere_pattern(mid+1, r)
            check = self._check_overlap(points, r_CG)
            check_high = self._check_overlap(points_high, r_CG)
            if check == False and check_high == True:
                points = _fast_sphere_pattern(mid, r)
                opt_points = mid
            elif check == True:
                max_points = mid - 1
            else:
                min_points = mid + 1

        for i, pos in enumerate(points):
            particle = mb.Compound(name="_CGN", pos=pos)
            self.add(particle, "_CGN[$]")

    def _check_overlap(self, points, radius):
        """ Determines if there is any overlap for a set of uniform spheres.

        Parameters:
        ----------
            points : np.ndarray (n, 3)
                Sphere locations
            radius : float
                Radius of spheres
        """
        dists = distance.cdist(points, points, 'euclidean')
        dists = dists[np.nonzero(dists)]

        return np.any(dists < 2.0 * radius)

class MME(mb.Compound):
    """ Coarse-grained alkane bead containing a CH2-CH2-CH3 group """
    def __init__(self):
        super(MME, self).__init__()
        self.add(mb.Particle(name="_MME"))
        self.add(mb.Port(anchor=self[0], orientation=[0, 1, 0], separation=0.15), 'up')

class MMM(mb.Compound):
    """ Coarse-grained alkane bead containing a CH2-CH2-CH2 group """
    def __init__(self):
        super(MMM, self).__init__()
        self.add(mb.Particle(name='_MMM'))
        self.add(mb.Port(anchor=self[0], orientation=[0, 1, 0], separation=0.15), 'up')
        self.add(mb.Port(anchor=self[0], orientation=[0, -1, 0], separation=0.15), 'down')

class CG_alkane(mb.Compound):
    """ A coarse-grained, alkane chain, with three carbons to a CG bead, optionally terminated by CH3 particles """
    def __init__(self, n=6, cap_front=True, cap_end=True):
        """ Initalize a CG_alkane Compound.

        Parameters:
        ----------
        n : int
            Number of backbone particles
        cap_front : boolean
            Add methyl gropu to beginning of chain ('down' port)
        cap_end : boolean
            Add methyl group to end of chain ('up' port)
        """

        super(CG_alkane, self).__init__()

        # Adjust length of Polymer for absense of methyl terminations.
        if not cap_front:
            n += 1
        if not cap_end:
            n += 1
        chain = mb.Polymer(MMM(), n=n-2, port_labels=('up', 'down'))
        self.add(chain, 'chain')
        if cap_front:
            self.add(MME(), 'methyl_front')
            mb.force_overlap(self['chain'], self['chain']['up'],
                             self['methyl_front']['up'])
        else:
            # Hoise part label to CG_alkane level.
            self.add(chain['up'], 'up', containment=False)

        if cap_end:
            self.add(MME(), 'methyl_end')
            mb.force_overlap(self['methyl_end'], self['methyl_end']['up'],
                             self['chain']['down'])
        else:
            # Hoist port label to CG_alkane level.
            self.add(chain['down'], 'down', containment=False)


class cgnp_patchy(mb.Compound):
    """
    Example class that would go in your recipe.

    Parameters
    ----------
    your_argument: int
        This is an example argument

    """
    def __init__(self):
        super(cgnp_patchy, self).__init__()
        nano = Nanoparticle(2.5, 0.6)
        self.add(nano, "nanoparticle")
