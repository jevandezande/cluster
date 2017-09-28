#!/usr/bin/env python3

import pytest
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
print(sys.path)
from cluster.cmolecule import *


class TestCMolecule:
    def test_init(self):
        cmol = CMolecule([['H', [0, 0, 0], 0],
                        ['O', [0, 0, 1],-1]])
        assert cmol.atoms == ['H', 'O']
        assert (cmol.xyz == [[0, 0, 0], [0, 0, 1]]).all()
        assert all(cmol.charges == [0, -1])
        assert cmol.charge == -1

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

        tmp_file = 'geom.xyz.tmp'
        cmol.write(tmp_file)
        cmol2 = CMolecule.read_from(tmp_file)
        assert cmol == cmol2
        os.remove(tmp_file)

    def test_other(self):
        pass
