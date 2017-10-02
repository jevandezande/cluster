#!/usr/bin/env python3
import os

import pytest

from pytest import approx
from cluster import Cluster
from cluster.cmolecule import CMolecule


class TestCluster:
    def test_init(self):
        cmol1 = CMolecule([['H', [0, 0, 0], 0],
                           ['N', [0, 0, 1], 1]])
        cmol2 = CMolecule([['H', [1, 0, 0], 0],
                           ['O', [1, 0, 1], -1]])
        cmol3 = CMolecule([['H', [2, 0, 0], 1],
                           ['F', [2, 0, 1], -2]])
        cluster = Cluster(cmol1, cmol2, cmol3)

        assert cluster.atoms == ['H', 'N', 'H', 'O', 'H', 'F']
        assert (cluster.xyz == [[0, 0, 0], [0, 0, 1], [1, 0, 0], [1, 0, 1], [2, 0, 0], [2, 0, 1]]).all()
        assert all(cluster.charges == [0, 1, 0, -1, 1, -2])
        assert cluster.charge == approx(-1)

        assert str(cluster) == """\
H       0.00000000    0.00000000    0.00000000  0.0000
N       0.00000000    0.00000000    1.00000000  1.0000
H       1.00000000    0.00000000    0.00000000  0.0000
O       1.00000000    0.00000000    1.00000000 -1.0000
H       2.00000000    0.00000000    0.00000000  1.0000
F       2.00000000    0.00000000    1.00000000 -2.0000"""

        assert len(cluster) == 6

        tmp_file = 'geom.xyz.tmp'
        cluster.write(tmp_file)
        cluster = Cluster.read_from(tmp_file, [2, 4])
        assert cluster == cluster
        os.remove(tmp_file)

    def test_write_input_file(self):
        cmol1 = CMolecule([['H', [0, 0, 0], 0],
                           ['N', [0, 0, 1], 1]])
        cmol2 = CMolecule([['H', [1, 0, 0], 0],
                           ['O', [1, 0, 1], -1]])
        cmol3 = CMolecule([['H', [2, 0, 0], 1],
                           ['F', [2, 0, 1], -2]])
        cluster = Cluster(cmol1, cmol2, cmol3)

        tmp_input_file = 'input.dat.tmp'

        options = {
            'ecp': 'SDD',
            'header': '! B3LYP def2-svp\n*xyz 0 1\n',
            'program': 'orca',
        }
        cluster.write_input_file(tmp_input_file, **options)

        with open(tmp_input_file) as f:
            inp_res = f.read()
        assert inp_res == """! B3LYP def2-svp
*xyz 0 1
    H               0.00000000    0.00000000    0.00000000
    N               0.00000000    0.00000000    1.00000000
    H>    0.0000    1.00000000    0.00000000    0.00000000 NewECP "SDD" end
    O>   -1.0000    1.00000000    0.00000000    1.00000000 NewECP "SDD" end
    Q     1.0000    2.00000000    0.00000000    0.00000000
    Q    -2.0000    2.00000000    0.00000000    1.00000000
*"""

        tmp_pc_file = 'input.pc.tmp'
        cluster.write_input_file(tmp_input_file, separate_pc=tmp_pc_file, **options)

        with open(tmp_input_file) as f:
            inp_res = f.read()
        assert inp_res == """! B3LYP def2-svp
*xyz 0 1
    H               0.00000000    0.00000000    0.00000000
    N               0.00000000    0.00000000    1.00000000
    H>    0.0000    1.00000000    0.00000000    0.00000000 NewECP "SDD" end
    O>   -1.0000    1.00000000    0.00000000    1.00000000 NewECP "SDD" end
*"""
        with open(tmp_pc_file) as f:
            pc_res = f.read()

        assert pc_res == """2
 1.0000    2.00000000    0.00000000    0.00000000
-2.0000    2.00000000    0.00000000    1.00000000
"""

        os.remove(tmp_input_file)
        os.remove(tmp_pc_file)

    def test_recharged(self):
        cmol1 = CMolecule([['H', [0, 0, 0], 0],
                           ['N', [0, 0, 1], 1]])
        cmol2 = CMolecule([['H', [1, 0, 0], 0],
                           ['O', [1, 0, 1], -1]])
        cmol3 = CMolecule([['H', [2, 0, 0], 1],
                           ['F', [2, 0, 1], -2]])
        cluster = Cluster(cmol1, cmol2, cmol3)

        recharged = cluster.recharged(3)

        assert recharged.charge == approx(3)
        assert recharged != cluster
        assert recharged.atoms == cluster.atoms
        assert (recharged.xyz == cluster.xyz).all()
        assert recharged.qc_mol == cluster.qc_mol
        assert recharged.br_mol.charges == approx(cluster.br_mol.charges + 1)
        assert recharged.pc_mol.charges == approx(cluster.pc_mol.charges + 1)

        recharged2 = cluster.recharged(2.9)
        assert recharged2 != recharged
        assert recharged2.charge == approx(2.9)
        assert recharged2.br_mol.charges == approx(cluster.br_mol.charges + 3.9/4)
        assert recharged2.pc_mol.charges == approx(cluster.pc_mol.charges + 3.9/4)
