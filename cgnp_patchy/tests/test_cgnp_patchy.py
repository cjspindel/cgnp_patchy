"""
Unit and regression test for the cgnp_patchy package.
"""

# Import package, test suite, and other packages as needed
import cgnp_patchy
import pytest
import sys
import mbuild as mb
import numpy as np

class BaseTest:
    @pytest.fixture
    def Core(self):
        return Nanoparticle(radius=2.5, bead_diameter=0.6) 
    
    @pytest.fixture
    def Alkane(self):
        return CGAlkane()

    @pytest.fixture
    def GraftedNanoparticle(self):
        return cgnp_patchy(radius=2.5, bead_diameter=0.6, chain_density=2.0)

class TestCGNPBuilder(BaseTest):
    def test_cgnp_patchy_imported():
        """ Sample test, will always pass so long as import statement worked """
        assert "cgnp_patchy" in sys.modules

    def test_import():
        """ Test that mBuild recipe import works """
        assert "cgnp_patchy" in vars(mb.recipes).keys()
