#!/usr/bin/env python3
import os

from cluster.molecule import Molecule
from cluster.crystal import Crystal, CubicCrystal, RectangularCrystal

import pytest
import numpy as np

from pytest import approx


class TestCrystal:
    def test_init(self):
        Crystal(1, 2, 3, 90, 90, 90, [['F', (0.5, 0.5, 0.5)]], 'oI')

        with pytest.raises(ValueError) as e:
            Crystal(1, 1, 1, 90, 90, 90, [['F', (0.5, 0.5, 0.5)]], 'AAA')

    def test_repr_str(self):
        geom = [['H', (0.5, 0.5, 0.5)], ['F', (0, 0, 0)]]
        cry = Crystal(1, 1, 1, 90, 90, 90, geom, 'cP')
        assert repr(cry) == '<Crystal(cP) HF>'
        geom_str = """cP, 1, 1, 1
H       0.50000000    0.50000000    0.50000000
F       0.00000000    0.00000000    0.00000000"""
        assert str(cry) == geom_str

    @pytest.mark.xfail(raises=NotImplementedError)
    def test_read_cif(self):
        Crystal.read_cif('abc.cif')

    def test_tile(self):
        geom = [['H', (0.5, 0.5, 0.5)], ['F', (0, 0, 0)]]
        cry = Crystal(1, 1, 1, 90, 90, 90, geom, 'cP')

        assert Molecule(geom) == cry.tile(0, 1, 0, 1, 0, 1)

        xyz = []
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    move = np.array([1, 1, 1]) * [i, j, k]
                    xyz.append([0.5, 0.5, 0.5] + move)
                    xyz.append([0.0, 0.0, 0.0] + move)
        mol = Molecule()
        mol.xyz = xyz
        mol.atoms = ['H', 'F'] * 8

        assert cry.tile(0, 2, 0, 2, 0, 2) == mol

    def test_tile_cuts(self):
        geom = [['H', (0.5, 0.5, 0.5)], ['F', (0, 0, 0)]]
        cry = Crystal(1, 1, 1, 90, 90, 90, geom, 'cP')

        # Cut at long distance (thus not cutting)
        cuts = [[(9, 9, 9), (1, 1, 1)]]
        assert Molecule(geom) == cry.tile(0, 1, 0, 1, 0, 1, cuts=cuts)

        # Cut multiple times at long distance (thus not cutting)
        cuts = [[(9, 9, 9), (1, 1, 1)], [(-9, -9, -9), (-1, -1, -1)]]
        assert Molecule(geom) == cry.tile(0, 1, 0, 1, 0, 1, cuts=cuts)

        # Cut everything
        cuts = [[(-9, -9, -9), (1, 1, 1)]]
        assert Molecule() == cry.tile(-1, 1, -1, 1, -1, 1, cuts=cuts)

        # Cut off extra tiling
        cuts = [[(1, 0, 0), (1, 0, 0)]]
        assert Molecule(geom) == cry.tile(0, 9, 0, 1, 0, 1, cuts=cuts)

        # Cut off plane
        cuts = [[(0, 1, 0), (0, 1, 0)]]
        assert len(cry.tile(0, 9, 0, 9, 0, 9, cuts=cuts)) == 162

        # Cut off plane2
        cuts = [[(0, 0, 0), (1, 1, 1)]]
        assert len(cry.tile(-2, 2, -2, 2, -2, 2, cuts=cuts)) == 76

        # Cut off all but the central atom
        cuts = [[(0.1, 0.1, 0.1), (1, 1, 1)], [(0.1, 0.1, -0.1), (1, 1, -1)], [(0.1, -0.1, 0.1), (1, -1, 1)], [(-0.1, 0.1, 0.1), (-1, 1, 1)], [(-0.1, -0.1, -0.1), (-1, -1, -1)]]
        tiled = cry.tile(-2, 2, -2, 2, -2, 2, cuts=cuts)
        assert Molecule([['F', (0, 0, 0)]]) == tiled


class TestCubicCrystal:
    def test_init(self):
        geom = [['H', (0.5, 0.5, 0.5)], ['F', (0, 0, 0)]]
        CubicCrystal(1, geom, 'cP')

    def test_repr(self):
        geom = [['H', (0.5, 0.5, 0.5)], ['F', (0, 0, 0)]]
        cry = CubicCrystal(1, geom, 'cP')
        assert repr(cry) == '<CubicCrystal(cP) HF>'


class RectangularCrystal:
    def test_init(self):
        geom = [['H', (0.5, 0.5, 0.5)], ['F', (0, 0, 0)]]
        RectangularCrystal(1, 2, 3, geom, 'cP')

    def test_repr(self):
        geom = [['H', (0.5, 0.5, 0.5)], ['F', (0, 0, 0)]]
        cry = RectangularCrystal(1, 2, 3, geom, 'cP')
        assert repr(cry) == '<RectangularCrystal(cP) HF>'
