import itertools
import numpy as np

from .molecule import Molecule
from collections import Counter


class Crystal:
    SUPPORTED_GROUPS = [
        'oP', 'oS', 'oI', 'oF',  # Orthorhombic
        'tP', 'tI',              # Tetragonal
        'cP', 'cI', 'cF',        # Cubic
    ]

    def __init__(self, a, b, c, alpha, beta, gamma, geom, space_group):
        """ A perfect crystal
        :params a, b, c, alpha, beta, gamma: crystal parameters
        :param geom: unique atoms and their locations within the cell
            Warning: all atoms must be within the region [0, a), [0, b), [0, c)
        :param space_group: the space group the crystal belongs to
            Warning: currently uses Pearson symbols, which do not uniquely
            identify the space group
            TODO: use Schönflies notation or some other method
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
        """ String of the geometry in xyz format """
        return f'{self.space_group}, {self.a}, {self.b}, {self.c}\n{self.molecule}'

    @property
    def atoms(self):
        """ Unique atoms in crystal """
        return self.molecule.atoms

    @property
    def xyz(self):
        """ Coordinates of unique atoms """
        return self.molecule.xyz

    @staticmethod
    def read_cif(file_name):
        """ Read a cif file """
        raise NotImplementedError()

        with open(file_name) as f:
            lines = f.readlines()

        for line in lines:
            if line[0] == '#':
                continue

    def tile(self, a1, a2, b1, b2, c1, c2, ratio=False, cut=None):
        """ Generates a molecule of the crystal with the given dimensions
        :params a1, a2, b1, b2, c1, c2: start and end in each direction
        :param ratio: convert from ratio to angstroms
        :param cut: two tuples of point and normal vector defining a plane
            WARNING: will eventually be changed to crystal face notation
        """
        num_a, num_b, num_c = a2 - a1, b2 - b1, c2 - c1

        xyz = self.xyz if not ratio else self.xyz*[self.a, self.b, self.c]
        tiled = Molecule()

        if self.space_group[0] in ['o', 't', 'c']:
            atoms = self.atoms*num_a*num_b*num_c
            tiled_xyz = np.zeros((len(atoms), 3))
            cell_size = np.array([self.a, self.b, self.c])
            position = 0
            # TODO: make vector dimensions num_a x num_b x num_c x 3 perform a single add
            for i in range(a1, a2):
                for j in range(b1, b2):
                    for k in range(c1, c2):
                        tiled_xyz[position:position + len(xyz), :] = xyz + cell_size*(i, j, k)
                        position += len(xyz)
        else:
            raise NotImplementedError(f'tile() is not implemented for {space_group}.')

        if cut is not None:
            (x0, y0, z0), (a, b, c) = cut
            cut_xyz = []
            cut_atoms = []
            for atom, (x, y, z) in zip(atoms, tiled_xyz):
                if a*(x - x0) + b*(y - y0) + c*(z - z0) < 0:
                    cut_xyz.append((x, y, z))
                    cut_atoms.append(atom)
            tiled_xyz = np.array(cut_xyz)
            atoms = cut_atoms

        tiled.xyz = tiled_xyz
        tiled.atoms = atoms

        return tiled


class CubicCrystal(Crystal):
    def __init__(self, a, geom, space_group):
        """ A cubic crystal (a=b=c)
        :param a: cube length
        :param geom: atoms and their locations within the cell
        :param space_group: the space group the crystal belongs to (cP, cI, cF)
            Warning: currently uses Pearson symbols, which do not uniquely
            identify the space group
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
        """ A crystal with all angles equal to 90°:
            orthorhombic (a≠b≠c) (oP, oS, oI, oF)
            tetragonal   (a=b≠c) (tP, tI)
            cubic        (a=b=c) (cP, cI, cF)
        The repeated structure is made purely by translations
        This mostly exists because I don't know much about space groups, and
        want a simplified system that works for me
        :params a, b, c: crystal parameters
        :param space_group: the space group the crystal belongs to
            Warning: currently uses Pearson symbols, which do not uniquely
            identify the space group
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
