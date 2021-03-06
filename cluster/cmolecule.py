import numpy as np

from .molecule import Molecule


class CMolecule(Molecule):
    def __init__(self, geom=None):
        """ Molecule with charges on all atoms
        :param geom: List of lists ordered as [[atom, [x, y, z], charge], ...]
        """
        self.geom = geom

    def __str__(self):
        """ A string of the geometry, filling out positions with zeros and
        spaces as needed
        """
        form = '{:<4}' + ' {:> 13.8f}' * 3 + ' {:>7.4f}'
        return '\n'.join([form.format(atom, *xyz, charge) for atom, xyz, charge in self])

    def __iter__(self):
        for atom, xyz, charge in zip(self.atoms, self.xyz, self.charges):
            yield atom, xyz, charge

    def __getitem__(self, i):
        """ Name, coordinates and charge of the ith atom """
        return (self.atoms[i], self.xyz[i], self.charges[i])

    def __setitem__(self, i, atom_xyz_charge):
        """ Sets the ith atom name, coordinates, and charge """
        atom, xyz, charge = atom_xyz_charge
        Molecule.check_atom(atom, xyz)
        self.atoms[i] = atom
        self.xyz[i] = xyz
        self.charges[i] = charge

    def __delitem__(self, i):
        """ Deletes the ith atom from the CMolecule """
        super().__delitem__(i)
        self.charges = np.delete(self.charges, i)

    def __eq__(self, other):
        """ Checks if CMolecules are (exactly) equivalent """
        return super().__eq__(other) and all(self.charges == other.charges)

    @staticmethod
    def from_Molecule(mol, charges):
        """ Generate a CMolecule from a Molecule
        :param mol: a Molecule
        :param charges: list of charges based on atom index or
            dictionary containing charges for every atom based on atom name
        """
        if isinstance(charges, list):
            geom = [(atom, xyz, charge) for (atom, xyz), charge in zip(mol, charges)]
        elif isinstance(charges, dict):
            geom = [(atom, xyz, charges[atom]) for atom, xyz in mol]
        else:
            raise ValueError(f'Expected list or dict for charges, got {type(charges)}')
        return CMolecule(geom)

    @property
    def charge(self):
        """ Overall charge of the molecule """
        return sum(self.charges)

    def insert(self, i, atom, xyz, charge):
        """ Insert the atom in the ith position """
        super().insert(i, atom, xyz)
        self.charges = np.insert(self.charges, i, charge)

    @property
    def geom(self):
        """ Return the geometry """
        return [[atom, list(xyz), charge] for atom, xyz, charge in self]

    @geom.setter
    def geom(self, geom):
        """ Set the geometry """
        self.atoms = []
        self.xyz = np.array([])
        self.charges = np.array([])

        if geom:
            atoms, xyzs, charges = zip(*geom)
            Molecule.check_geom(zip(atoms, xyzs))
            self.atoms = list(atoms)
            self.xyz = np.array(xyzs)
            self.charges = np.array(charges)

    def append(self, atom, xyz, charge):
        """ Append atom to CMolecule """
        super().append(atom, xyz)
        self.charges = np.append(self.charges, charge)

    @staticmethod
    def read_from(infile, charges=None):
        """ Read from a file
        :param infile: file to read from
        :param charges: list or dictionary of charges (else read from infile)
        """
        # Attempt to read as an XYZ file
        with open(infile) as f:
            lines = f.readlines()
        if lines[0].strip().isdigit():
            # Strip off length if provided
            lines = lines[2:]
        geom = []
        for i, line in enumerate(lines):
            if line.strip() == '':
                continue
            if charges is not None:
                atom, *xyz = line.split()[:4]
                charge = charges[atom] if isinstance(charges, dict) else charges[i]
            else:
                atom, *xyz, charge = line.split()[:5]
            x, y, z = map(float, xyz)
            geom.append([atom, np.array([x, y, z]), float(charge)])

        return CMolecule(geom)

    def write(self, outfile='geom.xyz', style='xyz', xyzlabel=True):
        """ Writes the geometry to the specified file
        :param outfile: file to write to
        :param style: xyz or latex style output
        :param xyzlabel: label xyz coordinates (i.e. put number of atoms at top)
        Prints the size at the beginning if desired (to conform to XYZ format)
        """
        out = ''
        if style == 'xyz':
            if xyzlabel:
                out += f'{len(self)}\n\n'
            out += str(self)
        elif style == 'latex':
            out += '\\begin{verbatim}\n' + str(self) + '\n\\end{verbatim}'
        else:
            raise SyntaxError('Invalid style')

        with open(outfile, 'w') as f:
            f.write(out)
