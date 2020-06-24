"""
Unit and regression test for the cgnp_patchy package.
"""

# Import package, test suite, and other packages as needed
import cgnp_patchy
import pytest
import sys
import mbuild as mb

class BaseTest:
    @pytest.fixture
    def Core(self):
        return Nanoparticle(radius=2.5, bead_diameter=0.6) 
    
    @pytest.fixture
    def Alkane(self):
        from cgnp_patchy.lib.chains import CGAlkane 
        return CGAlkane()

    @pytest.fixture
    def CGNanoparticle(self):
        from cgnp_patchy.cgnp_patchy import cgnp_patchy  
        np = cgnp_patchy(radius=2.5, bead_diameter=0.6, chain_density=2.0)
        return np 

class TestNanoparticleBuilder(BaseTest):
    def test_cgnp_patchy_imported(self):
        """ Sample test, will always pass so long as import statement worked """
        assert "cgnp_patchy" in sys.modules

    def test_import(self):
        """ Test that mBuild recipe import works """
        assert "cgnp_patchy" in vars(mb.recipes).keys()

    def test_chainlength(self):
        from cgnp_patchy.lib.chains import CGAlkane
        chain = CGAlkane(n=10)
        assert chain.n_particles == 10

    def test_save(self, CGNanoparticle):
        CGNanoparticle.save('nanoparticle.mol2', overwrite=True)

    def test_save_with_pattern(self):
        nanoparticle = mb.recipes.cgnp_patchy(radius=2.5, bead_diameter=0.6, chain_density=2.5, coating_pattern='bipolar')
        nanoparticle.save('bipolar_nanoparticle.mol2', overwrite=True)

    def test_save_with_backfill(self, Alkane):
        backfill_chain = Alkane
        nanoparticle = mb.recipes.cgnp_patchy(radius=2.5, bead_diameter=0.6, chain_density=2.5, coating_pattern='bipolar', backfill=backfill_chain)
        nanoparticle.save('backfill_nanoparticle.mol2', overwrite=True)
