"""
Unit and regression tests for coating patterns in the cgnp_patchy package.
"""

import pytest
import mbuild as mb
import numpy as np

class BaseTest():
    @pytest.fixture
    def IsotropicPattern(self):
        radius = 2.5
        point_density = 3.0
        pattern = mb.SpherePattern(int(point_density*4.0*np.pi*radius**2.0))
        return pattern

class TestPatterns(BaseTest):
    def test_mbuild_sphere_pattern(self, IsotropicPattern):
        radius = 2.5
        point_density = 3.0 
        pattern = mb.SpherePattern(int(point_density*4.0*np.pi*radius**2.0))
        assert len(pattern.points) == 235
        assert len(IsotropicPattern.points) == 235

    def test_bipolar_pattern(self):
        from cgnp_patchy.lib.patterns import BipolarPattern
        pattern = BipolarPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 189

    def test_equatorial_pattern(self):
        from cgnp_patchy.lib.patterns import EquatorialPattern
        pattern = EquatorialPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 188 

    def test_polar_pattern(self):
        from cgnp_patchy.lib.patterns import PolarPattern
        pattern = PolarPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 188 

    def test_cube_pattern(self):
        from cgnp_patchy.lib.patterns import CubePattern
        pattern = CubePattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 163

    def test_random_pattern(self):
        from cgnp_patchy.lib.patterns import RandomPattern
        pattern = RandomPattern(radius=2.5, chain_density=3.0, seed=123)
        assert len(pattern.points) > 200 and len(pattern.points) < 300

    def test_ring_pattern(self):
        from cgnp_patchy.lib.patterns import RingPattern
        pattern = RingPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 181 

    def test_square_pattern(self):
        from cgnp_patchy.lib.patterns import SquarePattern
        pattern = SquarePattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 185

    def test_tetrahedral_pattern(self):
        from cgnp_patchy.lib.patterns import TetrahedralPattern
        pattern = TetrahedralPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 163 
