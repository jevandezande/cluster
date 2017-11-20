from collections import Counter
from qgrep.molecule import Molecule

import itertools
import numpy as np


class Crystal:
    SUPPORTED_GROUPS = [
        'oP', 'oS', 'oI', 'oF',  # Orthorhombic
        'tP', 'tI',              # Tetragonal
        'cP', 'cI', 'cF',        # Cubic
    ]

    def __init__(self, a, b, c, alpha, beta, gamma, geom, space_group):
        """
        :params a, b, c, alpha, beta, gamma: crystal parameters
        :param geom: atoms and their locations within the cell
            Warning: all atoms must be within the region [0, a), [0, b), [0, c)
        :param space_group: the space group the crystal belongs
            Warning: currently uses Pearson symbols, which do not uniquely
            identify the space group
            TODO: use Schönflies notation
        """
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.molecule = Molecule(geom)

        counts = Counter(self.atoms)
        self.name = ''
        for atom, count in Counter(self.atoms).items():
            self.name += f'{atom}{count}' if count > 1 else f'{atom}'

        if space_group not in self.SUPPORTED_GROUPS:
            raise ValueError(f'space_group: {space_group} is invalid or not yet supported.')
        self.space_group = space_group

    def __repr__(self):
        return f'<Crystal({self.space_group}) {self.name}>'

    def __str__(self):
        """
        Return a string of the geometry in xyz format.
        """
        out = f'{self.space_group}, {self.a}, {self.b}, {self.c}\n'
        return out + str(self.molecule)

    @property
    def atoms(self):
        return self.molecule.atoms

    @property
    def xyz(self):
        return self.molecule.xyz

    @staticmethod
    def read_cif(file_name):
        raise NotImplementedError()

        with open(file_name) as f:
            lines = f.readlines()

        for line in lines:
            if line[0] == '#':
                continue

    def tile(self, a1, a2, b1, b2, c1, c2, ratio=False):
        """
        Generates a molecule of the crystal with the given dimensions
        :params num_a, num_b, num_c: number of times to repeat along the coordinate
        :param ratio: convert from ratio to angstroms
        """
        num_a, num_b, num_c = a2 - a1, b2 - b1, c2 - c1

        xyz = self.xyz if not ratio else self.xyz*[self.a, self.b, self.c]
        tiled = Molecule()

        if self.space_group[0] in ['o', 't', 'c']:
            tiled_xyz = np.zeros((len(xyz) * num_a * num_b * num_c, 3))
            position = 0
            vector = np.array([self.a, self.b, self.c])
            # TODO: make vector dimensions num_a x num_b x num_c x 3 perform a single add
            for i in range(a1, a2):
                for j in range(b1, b2):
                    for k in range(c1, c2):
                        tiled_xyz[position:position + len(xyz), :] = xyz + vector*(i, j, k)
                        position += len(xyz)

            tiled.xyz = tiled_xyz
            tiled.atoms = self.atoms*num_a*num_b*num_c
        else:
            raise NotImplementedError(f'tile() is not implemented for {space_group}.')

        return tiled


class CubicCrystal(Crystal):
    def __init__(self, a, geom, space_group):
        """
        A cubic crystal (a=b=c) (cP, cI, cF)
        :param a: cube length
        :param geom: atoms and their locations within the cell
        """
        if space_group not in ['cP', 'cI', 'cF']:
            raise ValueError('Invalid space group for CubicCrystal')
        super().__init__(a, a, a, 90, 90, 90, geom, space_group)

    def __repr__(self):
        return f'<CubicCrystal({self.space_group}) {self.name}>'


class RectangularCrystal(Crystal):
    SUPPORTED_GROUPS = [
        'oP', 'oS', 'oI', 'oF',  # Orthorhombic
        'tP', 'tI',              # Tetragonal
        'cP', 'cI', 'cF',        # Cubic
    ]

    def __init__(self, a, b, c, geom, space_group):
        """
        Any crystal with all angles equal to 90°:
            orthorhombic (a≠b≠c) (oP, oS, oI, oF)
            tetragonal   (a=b≠c) (tP, tI)
            cubic        (a=b=c) (cP, cI, cF)
        The repeated structure is made purely by translations
        This mostly exists because I don't know much about space groups, and
        want a simplified system that works for me
        :params a, b, c: crystal parameters
        """
        if space_group not in ['oP', 'oS', 'oI', 'oF',
                               'tP', 'tI',
                               'cP', 'cI', 'cF']:
            raise ValueError('Invalid space group for RectangularCrystal')
        super().__init__(a, b, c, 90, 90, 90, geom, space_group)

    def __repr__(self):
        return f'<RectangularCrystal({self.space_group}) {self.name}>'


if __name__ == '__main__':
    cry = Crystal(1, 1, 1, 90, 90, 90, [['F', [0.5, 0.5, 0.5]]], 'cP')
    cube = CubicCrystal(1, [['F', [0.5, 0.5, 0.5]]], 'cP')
    simple = RectangularCrystal(1, 1, 1, [['H', [0, 0, 0]], ['F', [0.5, 0.5, 0.5]]], 'cP')
    simple.tile(0, 2, 0, 2, 0, 2, ratio=True)
    geom = [
        ['Ir', [-2.25255,  0.00000,  0.78965]],
        ['Ir', [ 2.25255,  0.00000,  0.78965]],
        ['Ir', [ 0.00000, -2.25255, -0.78965]],
        ['Ir', [ 0.00000,  2.25255, -0.78965]],
        ['O' , [-0.86633,  1.38622,  0.78965]],
        ['O' , [ 0.86633, -1.38622,  0.78965]],
        ['O' , [ 0.86633,  3.11888,  0.78965]],
        ['O' , [-1.38622, -0.86633, -0.78965]],
        ['O' , [ 1.38622,  0.86633, -0.78965]],
        ['O' , [-3.11888,  0.86633, -0.78965]],
        ['O' , [ 3.11888, -0.86633, -0.78965]],
        ['O' , [ 3.63877,  1.38622,  0.78965]],
    ]
    IrO2 = RectangularCrystal(4.5051, 4.5051, 3.1586, geom, 'cP')
    IrO2.tile(0, 2, 0, 2, 0, 2)
