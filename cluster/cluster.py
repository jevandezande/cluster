import numpy as np
from .cmolecule import CMolecule


class Cluster:
    def __init__(self, qc_mol, br_mol, pc_mol):
        """ Cluster model for X-ray computations
        """
        self.qc_mol = qc_mol
        self.br_mol = br_mol
        self.pc_mol = pc_mol

    @staticmethod
    def generate_from_ranges(cmol, qc, br):
        assert qc >= 0
        assert qc <= br
        assert br <= len(cmol)
        return Cluster.generate_from_indices(cmol, range(qc), range(qc, qc + br))

    @staticmethod
    def generate_from_indices(cmol, qc, br):
        qc_geom = []
        br_geom = []
        pc_geom = []
        for i, (atom, xyz, charge) in enumerate(cmol):
            if i in qc:
                qc_geom.append([atom, xyz, float(charge)])
            elif i in br:
                br_geom.append([atom, xyz, float(charge)])
            else:
                pc_geom.append([atom, xyz, float(charge)])

        return Cluster(CMolecule(qc_geom), CMolecule(br_geom), CMolecule(pc_geom))

    def __str__(self):
        return '{}\n{}\n{}'.format(self.qc_mol, self.br_mol, self.pc_mol)

    def __len__(self):
        return len(self.qc_mol.atoms) + len(self.br_mol.atoms) + len(self.pc_mol.atoms)

    def __iter__(self):
        yield from self.qc_mol
        yield from self.br_mol
        yield from self.pc_mol

    @property
    def atoms(self):
        return self.qc_mol.atoms + self.br_mol.atoms + self.pc_mol.atoms

    @property
    def xyz(self):
        xyz = np.zeros((len(self), 3))

        xyz[:len(self.qc_mol), ...] = self.qc_mol.xyz
        xyz[len(self.qc_mol):len(self.qc_mol) + len(self.br_mol), ...] = self.br_mol.xyz
        xyz[-len(self.pc_mol):, ...] = self.pc_mol.xyz

        return xyz

    @property
    def charge(self):
        return sum(self.charges)

    @property
    def charges(self):
        charges = np.zeros(len(self))
        charges[:len(self.qc_mol), ...] = self.qc_mol.charges
        charges[len(self.qc_mol):len(self.qc_mol) + len(self.br_mol), ...] = self.br_mol.charges
        charges[-len(self.pc_mol):, ...] = self.pc_mol.charges

        return charges

    @staticmethod
    def read_from(infile, groups):
        """Read from a file
        :param infile: file to read from
        :param groups: Atom groupings corresponding to [qc, br]
            either an integer or a setlike object
            all other atoms placed in point charge
        """
        cmol = CMolecule.read_from(infile)
        if isinstance(groups[0], int):
           return Cluster.generate_from_ranges(cmol, *groups)
        return Cluster.generate_from_indices(cmol, *groups)

    def write(self, outfile, label=True, style='xyz'):
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
            out = '\\begin{verbatim}\n' + atoms + '\n\\end{verbatim}'
        else:
            raise SyntaxError('Invalid style')
        with open(outfile, 'w') as f:
            f.write(out)
