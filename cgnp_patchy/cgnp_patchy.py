"""
Primary function of recipe here
"""

import mbuild as mb

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
        # Sample of how a compound would be added in mBuild
        sample = mb.Particle(pos=[0.0, 0.0, 0.0], name='test')
        self.add(sample)
