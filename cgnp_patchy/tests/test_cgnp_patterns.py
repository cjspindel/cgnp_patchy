"""
Unit and regression tests for coating patterns in the cgnp_patchy package.
"""
import pytest
import mbuild as mb
import numpy as np
from cgnp_patchy.lib.utils.count_points import count_patch_points

class BaseTest():
    @pytest.fixture
    def IsotropicPattern(self):
        radius = 2.5
        point_density = 3.0
        pattern = mb.SpherePattern(int(point_density*4.0*np.pi*radius**2.0))
        pattern.scale(radius) 
        return pattern

class TestPatterns(BaseTest):
    ''' 
    Each pattern test checks if:
    - the number of points included in the pattern is equal to a set amount,
    - the number of points removed to create the pattern is equal to a set amount,
    - the number of points remaining and the number of points removed is equal to the number of points in 
    ''' 
    def test_mbuild_sphere_pattern(self, IsotropicPattern):
        ''' If this fails, it's likely mBuild's SpherePattern has been changed ''' 
        assert len(IsotropicPattern.points) == 235

    def test_bipolar_pattern(self, IsotropicPattern):
        from cgnp_patchy.lib.patterns import BipolarPattern
        pattern = BipolarPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 189
        assert count_patch_points(pattern, 2.5, 3.0) == 46
        assert len(pattern.points) + count_patch_points(pattern, 2.5, 3.0) == len(IsotropicPattern.points)

    def test_equatorial_pattern(self, IsotropicPattern):
        from cgnp_patchy.lib.patterns import EquatorialPattern
        pattern = EquatorialPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 188 
        assert count_patch_points(pattern, 2.5, 3.0) == 47 
        assert len(pattern.points) + count_patch_points(pattern, 2.5, 3.0) == len(IsotropicPattern.points)

    def test_polar_pattern(self, IsotropicPattern):
        from cgnp_patchy.lib.patterns import PolarPattern
        pattern = PolarPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 188 
        assert count_patch_points(pattern, 2.5, 3.0) == 47
        assert len(pattern.points) + count_patch_points(pattern, 2.5, 3.0) == len(IsotropicPattern.points)

    def test_cube_pattern(self, IsotropicPattern):
        from cgnp_patchy.lib.patterns import CubePattern
        pattern = CubePattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 163
        assert count_patch_points(pattern, 2.5, 3.0) == 72
        assert len(pattern.points) + count_patch_points(pattern, 2.5, 3.0) == len(IsotropicPattern.points)

    def test_random_pattern(self, IsotropicPattern):
        from cgnp_patchy.lib.patterns import RandomPattern
        pattern = RandomPattern(radius=2.5, chain_density=3.0, seed=123)
        assert len(pattern.points) > 200 and len(pattern.points) < 300

    def test_ring_pattern(self, IsotropicPattern):
        from cgnp_patchy.lib.patterns import RingPattern
        pattern = RingPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 181 
        assert count_patch_points(pattern, 2.5, 3.0) == 54
        assert len(pattern.points) + count_patch_points(pattern, 2.5, 3.0) == len(IsotropicPattern.points)

    def test_square_pattern(self, IsotropicPattern):
        from cgnp_patchy.lib.patterns import SquarePattern
        pattern = SquarePattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 185
        assert count_patch_points(pattern, 2.5, 3.0) == 50
        assert len(pattern.points) + count_patch_points(pattern, 2.5, 3.0) == len(IsotropicPattern.points)

    def test_tetrahedral_pattern(self, IsotropicPattern):
        from cgnp_patchy.lib.patterns import TetrahedralPattern
        pattern = TetrahedralPattern(radius=2.5, chain_density=3.0, fractional_sa=0.2)
        assert len(pattern.points) == 163
        assert count_patch_points(pattern, 2.5, 3.0) == 72
        assert len(pattern.points) + count_patch_points(pattern, 2.5, 3.0) == len(IsotropicPattern.points)

