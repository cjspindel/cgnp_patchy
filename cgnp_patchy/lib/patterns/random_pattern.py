from __future__ import division

import mbuild as mb
import numpy as np


class RandomPattern(mb.Pattern):
    """A pattern where points are distributed semi-randomly

    Parameters
    ----------
    chain_density : float
        Density of chain coating on the nanoparticle (chains / nm^2)
    radius : float
        Radius of the nanoparticle (nm)
    seed : int, optional, default=12345
        Seed for the random number generator
    """
    def __init__(self, chain_density, radius, seed=12345):
        np.random.seed(12345)
        pattern = mb.SpherePattern(int(chain_density * 20.0 * np.pi * radius**2.0))
        pattern.scale(radius)
        np.random.shuffle(pattern.points)
        points = pattern.points[:int(len(pattern.points)/5)]

        super(RandomPattern, self).__init__(points=points, orientations=None)
