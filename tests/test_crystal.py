#!/usr/bin/env python3
import os

import pytest

from pytest import approx

from cluster.crystal import Crystal, CubicCrystal, RectangularCrystal


class TestCrystal:
    def test_init(self):
        Crystal(1, 2, 3, 90, 90, 90, [['F', 0.5, 0.5, 0.5]], 'oI')

        with pytest.raises(ValueError) as e:
            Crystal(1, 1, 1, 90, 90, 90, [['F', 0.5, 0.5, 0.5]], 'AAA')

    def test_repr_str(self):
        geom = [['H', 0.5, 0.5, 0.5], ['F', 0, 0, 0]]
        cry = Crystal(1, 1, 1, 90, 90, 90, geom, 'cP')
        assert repr(cry) == '<Crystal(cP) HF>'
        geom_str = """\
H       0.50000000    0.50000000    0.50000000
F       0.00000000    0.00000000    0.00000000"""
        assert cry.geom_str() == geom_str

    @pytest.mark.xfail(raises=NotImplementedError)
    def test_read_cif(self):
        Crystal.read_cif('abc.cif')

    def test_tile(self):
        geom = [['H', 0.5, 0.5, 0.5], ['F', 0, 0, 0]]
        cry = Crystal(1, 1, 1, 90, 90, 90, geom, 'cP')
        assert cry.geom_str() == cry.tile(0, 1, 0, 1, 0, 1)
        tile_str_2_2_2 = """\
H       0.50000000    0.50000000    0.50000000
F       0.00000000    0.00000000    0.00000000
H       0.50000000    0.50000000    1.50000000
F       0.00000000    0.00000000    1.00000000
H       0.50000000    1.50000000    0.50000000
F       0.00000000    1.00000000    0.00000000
H       0.50000000    1.50000000    1.50000000
F       0.00000000    1.00000000    1.00000000
H       1.50000000    0.50000000    0.50000000
F       1.00000000    0.00000000    0.00000000
H       1.50000000    0.50000000    1.50000000
F       1.00000000    0.00000000    1.00000000
H       1.50000000    1.50000000    0.50000000
F       1.00000000    1.00000000    0.00000000
H       1.50000000    1.50000000    1.50000000
F       1.00000000    1.00000000    1.00000000"""
        assert cry.tile(0, 2, 0, 2, 0, 2) == tile_str_2_2_2


class TestCubicCrystal:
    def test_init(self):
        geom = [['H', 0.5, 0.5, 0.5], ['F', 0, 0, 0]]
        CubicCrystal(1, geom, 'cP')

    def test_repr(self):
        geom = [['H', 0.5, 0.5, 0.5], ['F', 0, 0, 0]]
        cry = CubicCrystal(1, geom, 'cP')
        assert repr(cry) == '<CubicCrystal(cP) HF>'


class RectangularCrystal:
    def test_init(self):
        geom = [['H', 0.5, 0.5, 0.5], ['F', 0, 0, 0]]
        RectangularCrystal(1, 2, 3, geom, 'cP')

    def test_repr(self):
        geom = [['H', 0.5, 0.5, 0.5], ['F', 0, 0, 0]]
        cry = RectangularCrystal(1, 2, 3, geom, 'cP')
        assert repr(cry) == '<RectangularCrystal(cP) HF>'
