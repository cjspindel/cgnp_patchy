import mbuild as mb
import numpy as np

from cgnp_patchy.lib.nanoparticles import Nanoparticle
from cgnp_patchy.lib.chains import CG_alkane
from cgnp_patchy.lib.patterns import *

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
    coating_pattern : str, optional, default='isotropic'
        Type of pattern for the chain coating.
        Supported types are 'polar', 'bipolar', 'isotropic', 'equatorial', 'square', 'random', 'cube', 'tetrahedral', and 'ring'.
    fractional_sa : float, default=0.2
        Fractional surface rea of the nanoparticle to exclude coating (nm^2)
    """
    def __init__(self, radius, bead_diameter, chain_density, backfill=None, coating_pattern='isotropic', fractional_sa=0.2, **kwargs):
        super(cgnp_patchy, self).__init__()

        self.bead_diameter = bead_diameter
        
        chain = CG_alkane()
        nano = Nanoparticle(radius, bead_diameter)
        self.add(nano, 'nanoparticle')

        isotropic_pattern = mb.SpherePattern(int(chain_density*4.0*np.pi*radius**2.0))
        isotropic_pattern.scale(radius)
       
        if coating_pattern == 'isotropic':
            pattern = isotropic_pattern
        elif coating_pattern == 'polar':
            pattern = PolarPattern(chain_density, radius, fractional_sa, **kwargs)
        elif coating_pattern == 'bipolar':
            pattern = BipolarPattern(chain_density, radius, fractional_sa, **kwargs)
        elif coating_pattern == 'equatorial':
            pattern = EquatorialPattern(chain_density, radius, fractional_sa, **kwargs)
        elif coating_pattern == 'square':
            pattern = SquarePattern(chain_density, radius, fractional_sa, **kwargs)
        elif coating_pattern == 'random':
            pattern = RandomPattern(chain_density, radius, fractional_sa, **kwargs)
        elif coating_pattern == 'cube':
            pattern = CubePattern(chain_density, radius, fractional_sa, **kwargs)
        elif coating_pattern == 'tetrahedral':
            pattern = TetrahedralPattern(chain_density, radius, fractional_sa, **kwargs)
        elif coating_pattern == 'ring':
            pattern = RingPattern(chain_density, radius, fractional_sa, **kwargs)
        else:
            raise Exception("Coating pattern '{}' not supported. Valid options are 'polar', 'bipolar', 'isotropic', 'equatorial', 'square', 'random', 'cube', 'tetrahedral', and 'ring'.".format(coating_pattern))

        if backfill and coating_pattern == 'random':
            raise Exception("Backfill not supported for coating pattern type 'random'.")

        if backfill:
            backfill_points = []
            for point in isotropic_pattern.points:
                if not np.all(np.isin(point, isotropic_pattern.points)):
                    backfill_points.append(point)
        
        # Hacky workaround until apply_to_compound below can be used
        for pos in pattern.points:
            port = mb.Port(anchor=self['nanoparticle'], orientation=pos, separation=radius)
            self['nanoparticle'].add(port, 'port[$]')
            self.add(chain)
            mb.force_overlap(chain, chain['up'], port)  

        #chain_protos, empty_backfill = isotropic_pattern.apply_to_compound(guest=chain, guest_port_name='up', host=self['nanoparticle'])
        #self.add(chain_protos)

        
        if backfill:
            pattern.points = np.array(backfill_points)
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
