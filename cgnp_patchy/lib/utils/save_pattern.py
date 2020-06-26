import mbuild as mb
import numpy as np

def save_pattern(filename, pattern, overwrite=False):
    lj_proto = mb.Compound(name='LJ')
    lj_box = mb.Compound()
    for pos in pattern:
        lj_particle = mb.clone(lj_proto)
        lj_particle.translate(pos)
        lj_box.add(lj_particle)
    lj_box.save(filename, overwrite=overwrite)
