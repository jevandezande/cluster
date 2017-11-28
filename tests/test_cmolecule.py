#!/usr/bin/env python3
import os
import sys

import pytest

from pytest import approx

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from cluster.cmolecule import *


class TestCMolecule:
    def test_init(self):
        geom = [['H', (0, 0, 0), 0],
                ['O', (0, 0, 1), -1]]
        cmol = CMolecule(geom)
        assert cmol.atoms == ['H', 'O']
        assert (cmol.xyz == [[0, 0, 0], [0, 0, 1]]).all()
        assert all(cmol.charges == [0, -1])
        assert cmol.charge == approx(-1)

        assert str(cmol) == """\
H       0.00000000    0.00000000    0.00000000  0.0000
O       0.00000000    0.00000000    1.00000000 -1.0000"""

        assert len(cmol) == 2

        cmol[0] = 'B', [0, 1, 2], 99
        assert cmol[0][0] == 'B'

        del cmol[0]
        assert len(cmol) == 1

        cmol.append('W', [4, 5, 6], -8)
        assert len(cmol) == 2
        assert cmol.charge == -9

    def test_read_from(self):
        geom = [['H', (0, 0, 0), 0],
                ['O', (0, 0, 1), -1]]
        cmol = CMolecule(geom)
        tmp_file = 'geom.xyz.tmp'
        cmol.write(tmp_file)
        cmol2 = CMolecule.read_from(tmp_file)
        assert cmol == cmol2
        cmol3 = CMolecule.read_from(tmp_file, charges=[5, -8])
        assert cmol != cmol3
        assert cmol3.charges == approx([5, -8])

        cmol4 = CMolecule.read_from(tmp_file, charges={'H': 5, 'O': -8})
        assert cmol4.charges == approx([5, -8])
        assert cmol3 == cmol4

        os.remove(tmp_file)

    def test_from_Molecule(self):
        geom = [['H', (0, 0, 0)],
                ['O', (0, 0, 1)]]
        mol = Molecule(geom)
        cmol = CMolecule.from_Molecule(mol, [0, -1])
        assert cmol.atoms == mol.atoms
        assert approx(cmol.xyz) == mol.xyz
        assert approx(cmol.charges) == [0, -1]
        cmol = CMolecule.from_Molecule(mol, {'H': 0, 'O': -1})
        assert cmol.atoms == mol.atoms
        assert approx(cmol.charges) == [0, -1]
        with pytest.raises(ValueError) as e:
            CMolecule.from_Molecule(mol, 'HO')

