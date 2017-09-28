#!/usr/bin/env python3
import os

import pytest

from cluster.cluster import Cluster
from cluster.cmolecule import CMolecule


class TestCluster:
    def test_init(self):
        cmol1 = CMolecule([['H', [0, 0, 0], 0],
                           ['N', [0, 0, 1], 1]])
        cmol2 = CMolecule([['H', [1, 0, 0], 0],
                           ['O', [1, 0, 1],-1]])
        cmol3 = CMolecule([['H', [2, 0, 0], 0],
                           ['F', [2, 0, 1], 0]])
        cluster = Cluster(cmol1, cmol2, cmol3)
        assert cluster.atoms == ['H', 'N', 'H', 'O', 'H', 'F']
        assert (cluster.xyz == [[0, 0, 0], [0, 0, 1], [1, 0, 0], [1, 0, 1], [2, 0, 0], [2, 0, 1]]).all()
        assert all(cluster.charges == [0, 1, 0, -1, 0, 0])
        assert cluster.charge == 0

        assert str(cluster) == """\
H       0.00000000    0.00000000    0.00000000  0.0000
N       0.00000000    0.00000000    1.00000000  1.0000
H       1.00000000    0.00000000    0.00000000  0.0000
O       1.00000000    0.00000000    1.00000000 -1.0000
H       2.00000000    0.00000000    0.00000000  0.0000
F       2.00000000    0.00000000    1.00000000  0.0000"""

        assert len(cluster) == 6

        tmp_file = 'geom.xyz.tmp'
        cluster.write(tmp_file)
        cluster = Cluster.read_from(tmp_file, [2, 4])
        assert cluster == cluster
        os.remove(tmp_file)

    def test_other(self):
        pass
