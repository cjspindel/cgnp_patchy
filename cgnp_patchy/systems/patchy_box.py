from __future__ import division

import os
from pkg_resources import resource_filename

import mbuild as mb
import numpy as np
import random
from copy import deepcopy

class PatchyBox(mb.Compound):
    def __init__(self, nano, n, box, seed=12345):
        super(PatchyBox, self).__init__()
        
        if type(nano) is not list:
            nano = [nano]
        if type(n) is not list:
            n = [n] 
        
        # Define positions for nanoparticles (use points to speed
        # this up)
        point_rep = mb.Particle(name='point0')
        d_vdw = max([max(nano[0].xyz[:,dim]) - min(nano[0].xyz[:,dim])
                     for dim in range(3)]) 
        point_box = mb.fill_box(point_rep, n[0], box, overlap=d_vdw+0.5, edge=d_vdw/2+0.25, seed=seed)
        
        for i, (np_proto, np_n) in enumerate(zip(nano[1:], n[1:])):
            point_rep = mb.Particle(name='point{:d}'.format(i+1))
            d_vdw = max([max(np_proto.xyz[:,dim]) - min(np_proto.xyz[:,dim])
                     for dim in range(3)])
            point_box = mb.solvate(solute=point_box, solvent=point_rep, n_solvent=np_n, box=box,
                                   overlap=d_vdw+0.5, edge=d_vdw/2+0.25, seed=seed)
        self.periodicity = box.lengths
        random.seed(seed)
        
        # Replicate the nanoparticle at the defined positions
        for particle in point_box.particles():
            nano_index = int(particle.name.strip('point'))
            nano_clone = mb.clone(nano[nano_index])
            nano_clone.spin(random.random()*2*np.pi, np.array([random.random(), random.random(), random.random()]) - 0.5)
            #mb.translate_to(nano_clone, particle.pos)
            nano_clone.translate_to(particle.pos)
            self.add(nano_clone)

if __name__ == "__main__":
    import mbuild as mb
    from cgnp_patchy.cgnp_patchy import cgnp_patchy
    from cgnp_patchy.lib.chains import CGAlkane
    chain_proto = CGAlkane(n=6, cap_front=False, cap_end=True)
    tnp = cgnp_patchy(radius=2.5, chain=chain_proto, chain_density=2.5, coating_pattern='bipolar', fractional_sa = 0.2)
    box = mb.Box(lengths=np.ones(3)*20)
    patchy_box = PatchyBox(tnp, n=10, box=box)
    
    forcefield_dir = resource_filename('cgnp_patchy', '/forcefields')
    patchy_box.save('patchy-box.hoomdxml',
                    forcefield_files=[os.path.join(forcefield_dir, 'cg-alkane.xml'),
                                      os.path.join(forcefield_dir, 'nano-0.6.xml')],
                    ref_distance=3.95, ref_energy=0.091493, overwrite=True)
