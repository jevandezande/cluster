import numpy as np
from qgrep.molecule import Molecule


class CMolecule(Molecule):
    def __init__(self, geom=None):
        """ Simple molecule class

        :param geom: List of lists ordered as [[atom, [x, y, z], charge], ...]
        """
        self.geom = geom

    def __str__(self):
        """
        Returns a string of the geometry, filling out positions with zeros and
        spaces as needed
        """
        form = '{:<4}' + ' {:> 13.8f}' * 3 + ' {:>7.4f}'
        return '\n'.join([form.format(atom, *xyz, charge) for atom, xyz, charge in self])

    def __iter__(self):
        for atom, xyz, charge in zip(self.atoms, self.xyz, self.charges):
            yield atom, xyz, charge

    def __getitem__(self, i):
        """Returns the ith atom"""
        return (self.atoms[i], self.xyz[i], self.charges[i])

    def __setitem__(self, i, atom_xyz_charge):
        """Sets the ith atom"""
        atom, xyz, charge = atom_xyz_charge
        Molecule.check_atom(atom, xyz)
        self.atoms[i] = atom
        self.xyz[i] = xyz
        self.charges[i] = charge

    def __delitem__(self, i):
        """Deletes the ith atom"""
        super().__delitem__(i)
        self.charges = np.delete(self.charges, i)

    def __eq__(self, other):
        return super().__eq__(other) and all(self.charges == other.charges)

    @property
    def charge(self):
        return sum(self.charges)

    def insert(self, i, atom, xyz, charge):
        """Insert the atom in the specified position"""
        super().insert(i, atom, xyz)
        self.charges = np.insert(self.charges, i, charge)

    @property
    def geom(self):
        """Return the geometry"""
        return [[atom, list(xyz), charge] for atom, xyz, charge in self]

    @geom.setter
    def geom(self, geom):
        """Set the geometry"""
        self.charges = np.array([])

        if geom:
            atoms, xyzs, charges = zip(*geom)
            Molecule.check_geom(zip(atoms, xyzs))
            self.atoms = list(atoms)
            self.xyz = np.array(xyzs)
            self.charges = np.array(charges)

    def append(self, atom, xyz, charge):
        """Append atom to geometry"""
        super().append(atom, xyz)
        self.charges = np.append(self.charges, charge)

    @staticmethod
    def read_from(infile):
        """Read from a file"""
        # Attempt to read as an XYZ file
        with open(infile) as f:
            lines = f.readlines()
        if lines[0].strip().isdigit():
            # Strip off length if provided
            lines = lines[2:]
        geom = []
        for line in lines:
            if line.strip() == '':
                continue
            atom, x, y, z, charge = line.split()[:5]
            geom.append([atom, np.array([float(x), float(y), float(z)]), float(charge)])

        return CMolecule(geom)

    def write(self, outfile="geom.xyz", label=True, style='xyz'):
        """
        Writes the geometry to the specified file
        Prints the size at the beginning if desired (to conform to XYZ format)
        """
        out = ''
        if style == 'xyz':
            if label:
                out += '{}\n\n'.format(len(self))
            out += str(self)
        elif style == 'latex':
            header = '{}\\\\\n'.format(len(self))
            line_form = '{:<2}' + ' {:> 13.6f}' * 3 + ' {:>7.4f}'
            atoms = [line_form.format(atom, *pos, charge) for atom, xyz, charge in self]
            atoms = '\n'.join(atoms)
            #out = header + '\\begin{verbatim}\n' + atoms + '\n\\end{verbatim}'
            out = '\\begin{verbatim}\n' + atoms + '\n\\end{verbatim}'
        else:
            raise SyntaxError('Invalid style')
        with open(outfile, 'w') as f:
            f.write(out)
