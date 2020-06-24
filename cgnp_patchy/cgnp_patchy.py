import mbuild as mb
import numpy as np

from cgnp_patchy.lib.nanoparticles import Core
from cgnp_patchy.lib.chains import CGAlkane
from cgnp_patchy.lib.patterns import *

class cgnp_patchy(mb.Compound):
    """
    Builds a tethered, coarse-grained nanoparticle.

    Parameters
    ----------
    radius : float
        Radius of the nanoparticle (nm)
    chain : mb.Compound
        Prototype of alkane chain to attach to the nanoparticle core. Not currently implemented.
    chain_density : float
        Density of chain coating on the nanoparticle (chains / nm^2)
    bead_diameter : float, default=0.6
        Diameter of CG particles in the nanoparticle core (nm)
    backfill : mb.Compound, optional, default=None
        Protoype of backfill to place at vacant sites on the nanoparticle
    coating_pattern : str, optional, default='isotropic'
        Type of pattern for the chain coating.
        Supported types are 'polar', 'bipolar', 'isotropic', 'equatorial', 'square', 'random', 'cube', 'tetrahedral', and 'ring'.
    fractional_sa : float, default=0.2
        Fractional surface rea of the nanoparticle to exclude coating (nm^2)
    """
    def __init__(self, radius, chain_density, bead_diameter=0.6, backfill=None, coating_pattern='isotropic', fractional_sa=0.2, **kwargs):
        super(cgnp_patchy, self).__init__()
        
        self.bead_diameter = bead_diameter
        
        nano = Core(radius, bead_diameter)
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
            pattern = RandomPattern(chain_density, radius, **kwargs)
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
        elif backfill and coating_pattern == 'isotropic':
            raise Exception("Backfill not supported for coating pattern type 'isotropic'.")

        if backfill:
            backfill_points = []
            for point in isotropic_pattern.points:
                if not np.all(np.isin(point, pattern.points)):
                    backfill_points.append(point)
        
        # Hacky workaround until apply_to_compound below can be used
        for pos in pattern.points:
            port = mb.Port(anchor=self['nanoparticle'], orientation=pos, separation=radius)
            self['nanoparticle'].add(port, 'port[$]')
            chain = CGAlkane()
            self.add(chain)
            mb.force_overlap(chain, chain['up'], port)  
        #chain_protos, empty_backfill = isotropic_pattern.apply_to_compound(guest=chain, guest_port_name='up', host=self['nanoparticle'])
        #self.add(chain_protos)
       
        if backfill:
            pattern.points = np.array(backfill_points)
            # Problems with apply_to_compound again, temporarily replaced with workaround used with pattern
            for pos in pattern.points:
                port = mb.Port(anchor=self['nanoparticle'], orientation=pos, separation=radius)
                self['nanoparticle'].add(port, 'b_port[$]')
                b_chain = mb.clone(backfill)
                self.add(b_chain)
                mb.force_overlap(b_chain, b_chain['up'], port)
            #backfill_protos, empty_backfill = pattern.apply_to_compound(backfill, guest_port_name='up', host=self['nanoparticle'])
            #self.add(backfill_protos)
       

        self.label_rigid_bodies(rigid_particles='_CGN')

        # This is a temporary workaround until the 'apply_to_compound' method in mBuild is fixed -Andrew
        # Has the problem this was working around been fixed yet? If so, this code can be updated.
        for bond in self.bonds():
            if bond[0].name == 'Core' or bond[1].name == 'Core':
                if bond[0].name == 'Core':
                    bond[1].rigid_id = 0
                if bond[1].name == 'Core':
                    bond[0].rigid_id = 0
                self.remove_bond(bond)
