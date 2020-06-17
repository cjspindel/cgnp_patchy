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


class cgnp_patchy(mb.Compound):
    """
    Builds a tethered, coarse-grained nanoparticle.

    Parameters
    ----------
    radius : float
        Radius of the nanoparticle (nm)
    bead_diameter : float
        Diameter of CG particles in the nanoparticle core (nm)
    chain : mb.Compound
        Prototype of alkane chain to attach to the nanoparticle core
    chain_density : float
        Density of chain coating on the nanoparticle (chains / nm^2)
    backfill : mb.Compound, optional, default=None
        Protoype of backfill to place at vacant sites on the nanoparticle

    """
    def __init__(self, radius, bead_diameter, chain_density, backfill=None, coating_pattern='isotropic', **kwargs):
        super(cgnp_patchy, self).__init__()

        self.bead_diameter = bead_diameter
        chain = CG_alkane()

        nano = Nanoparticle(radius, bead_diameter)
        self.add(nano, "nanoparticle")

        isotropic_pattern = mb.SpherePattern(int(chain_density*4.0*np.pi*radius**2.0))
        isotropic_pattern.scale(radius)

        if backfill:
            backfill_points = []
            for point in isotropic_pattern.points:
                if not np.all(np.isin(point, isotropic_pattern.points)):
                    backfill_points.append(point)
        for pos in isotropic_pattern.points:
            port = mb.Port(anchor=self['nanoparticle'], orientation=pos, separation=radius)
            self['nanoparticle'].add(port, "port[$]")

        chain_protos, empty_backfill = isotropic_pattern.apply_to_compound(guest=chain, guest_port_name='up', host=self['nanoparticle'])
        self.add(chain_protos)

        if backfill:
            isotropic_pattern.points = np.array(backfill_points)
            for pos in pattern.points:
                port = mb.Port(anchor=self['nanoparticle'], orientation=pos, separation=radius)
            backfill_protos, empty_backfill = pattern.apply_to_compound(backfill, guest_port_name='up', host=self['nanoparticle'])
            self.add(backfill_protos)

        self.label_rigid_bodies(rigid_particles='_CGN')

        # This is a temporary workaround until the 'apply_to_compound' method in mBuild is fixed -Andrew
        # Has the problem this was working around been fixed yet? If so, this code can be updated.
        for bond in self.bonds():
            if bond[0].name == 'Nanoparticle' or bond[1].name == 'Nanoparticle':
                if bond[0].name == 'Nanoparticle':
                    bond[1].rigid_id = 0
                if bond[1].name == 'Nanoparticle':
                    bond[0].rigid_id = 0
                self.remove_bond(bond)
