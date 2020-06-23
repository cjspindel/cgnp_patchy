"""
Unit and regression test for the cgnp_patchy package.
"""

# Import package, test suite, and other packages as needed
import cgnp_patchy
import pytest
import sys
import mbuild as mb
import numpy as np
from cgnp_patchy.lib.chains import CGAlkane

'''
class BaseTest:
    @pytest.fixture
    def Core(self):
        return Nanoparticle(radius=2.5, bead_diameter=0.6) 
    
    @pytest.fixture
    def Alkane(self):
        return CGAlkane()

    @pytest.fixture
    def CGNanoparticle(self):
        from cgnp_patchy.cgnp_patchy import cgnp_patchy  
        np = cgnp_patchy(radius=2.5, bead_diameter=0.6, chain_density=2.0)
        return np 
'''

def test_cgnp_patchy_imported():
    """ Sample test, will always pass so long as import statement worked """
    assert "cgnp_patchy" in sys.modules

def test_import():
    """ Test that mBuild recipe import works """
    assert "cgnp_patchy" in vars(mb.recipes).keys()

def test_alkane():
    chain = CGAlkane() 
    return chain

def test_chainlength():
    chain = CGAlkane(n=10)
    assert chain.n_particles == 10

def test_save():
    nanoparticle = mb.recipes.cgnp_patchy(radius=2.5, bead_diameter=0.6, chain_density=2.0) 
    nanoparticle.save('nanoparticle.mol2', overwrite=True)

