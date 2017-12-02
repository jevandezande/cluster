import numpy as np
import itertools
from .cmolecule import CMolecule


class Cluster:
    def __init__(self, qc_mol, br_mol, pc_mol):
        """ Cluster model for X-ray computations
        All mols must be CMolecules
        :param qc_mol: quantum region
        :param br_mol: boundary region (where a capped ECP will be placed)
        :param pc_mol: point charge region
        """
        assert isinstance(qc_mol, CMolecule)
        assert isinstance(br_mol, CMolecule)
        assert isinstance(pc_mol, CMolecule)

        self.qc_mol = qc_mol
        self.br_mol = br_mol
        self.pc_mol = pc_mol

    @staticmethod
    def generate_from_ranges(cmol, qc, br):
        """ Generate a cluster from a CMolecule and starting indices
        :params qc, br: start of the corresponding regions
        """
        assert qc >= 0
        assert qc <= br
        assert br <= len(cmol)
        return Cluster.generate_from_indices(cmol, range(qc), range(qc, qc + br))

    @staticmethod
    def generate_from_indices(cmol, qc, br):
        """ Generate a cluster from a CMolecule and indices
        :params qc, br: corresponding regions
        Everything else is point charges.
        """
        assert all([q <= len(cmol) for q in qc])
        assert all([q <= len(cmol) for q in br])

        qc_geom = []
        br_geom = []
        pc_geom = []
        overlap = set(qc) & set(br)
        if overlap:
            raise ValueError(f'Overlapping regions: {overlap}')

        for i, (atom, xyz, charge) in enumerate(cmol):
            if i in qc:
                qc_geom.append([atom, xyz, float(charge)])
            elif i in br:
                br_geom.append([atom, xyz, float(charge)])
            else:
                pc_geom.append([atom, xyz, float(charge)])

        return Cluster(CMolecule(qc_geom), CMolecule(br_geom), CMolecule(pc_geom))

    def __repr__(self):
        """ Representation of cluster containing the number of atoms in each region """
        return f'<Cluster {len(self.qc_mol)}/{len(self.br_mol)}/{len(self.pc_mol)}>'

    def __str__(self):
        """ String of the geometry in xyz style, filling out positions with
        zeros and spaces as needed
        """
        return self.str(sep=False)

    def str(self, sep=False):
        """ String of the geometry in xyz style, filling out positions with zeros
        and spaces as needed
        :param sep: add a separator between sections
        """
        dashes = '-'*54 + '\n' if sep else ''
        return f'{self.qc_mol}\n{dashes}{self.br_mol}\n{dashes}{self.pc_mol}'

    def __len__(self):
        """
        Total number of atoms
        """
        return len(self.qc_mol.atoms) + len(self.br_mol.atoms) + len(self.pc_mol.atoms)

    def __iter__(self):
        """ Iterate over all atoms """
        yield from self.qc_mol
        yield from self.br_mol
        yield from self.pc_mol

    @property
    def atoms(self):
        """ List of all atoms """
        return self.qc_mol.atoms + self.br_mol.atoms + self.pc_mol.atoms

    @property
    def xyz(self):
        """ np.array of all coordinates """
        xyz = np.zeros((len(self), 3))

        xyz[:len(self.qc_mol), ...] = self.qc_mol.xyz
        xyz[len(self.qc_mol):len(self.qc_mol) + len(self.br_mol), ...] = self.br_mol.xyz
        xyz[-len(self.pc_mol):, ...] = self.pc_mol.xyz

        return xyz

    @property
    def charge(self):
        """ Total charge of the system """
        return sum(self.charges)

    @property
    def charges(self):
        """ np.array of the charges on the atoms
        Note: this will include the charge on the qc region atoms as well
        """
        charges = np.zeros(len(self))
        charges[:len(self.qc_mol), ...] = self.qc_mol.charges
        charges[len(self.qc_mol):len(self.qc_mol) + len(self.br_mol), ...] = self.br_mol.charges
        charges[-len(self.pc_mol):, ...] = self.pc_mol.charges

        return charges

    @staticmethod
    def read_from(infile, groups, charges=None):
        """ Read from a file
        :param infile: file to read from
        :param groups: Atom groupings corresponding to [qc, br]
            either an integer or a setlike object
            all other atoms placed in point charge
        :param charges: list or dictionary of charges (else read from infile)
        """
        cmol = CMolecule.read_from(infile, charges)
        if isinstance(groups[0], int):
            return Cluster.generate_from_ranges(cmol, *groups)
        return Cluster.generate_from_indices(cmol, *groups)

    def write(self, outfile, style='xyz', xyzlabel=True):
        """ Write to a file
        :param outfile: file to write to
        :param style: xyz or latex style output
        :param xyzlabel: label xyz coordinates (i.e. put number of atoms at top)
        """
        out = ''
        if style == 'xyz':
            if xyzlabel:
                out += f'{len(self)}\n\n'
            out += str(self)
        elif style == 'latex':
            header = f'{len(self)}\\\\\n'
            line_form = '{:<2}' + ' {:> 13.6f}' * 3 + ' {:>7.4f}'
            atoms = [line_form.format(atom, *pos, charge) for atom, xyz, charge in self]
            atoms = '\n'.join(atoms)
            out = '\\begin{verbatim}\n' + atoms + '\n\\end{verbatim}'
        else:
            raise SyntaxError('Invalid style')
        with open(outfile, 'w') as f:
            f.write(out)

    def write_input_file(self, outfile='input.dat', **options):
        """ Write an input file
        Options
        program: what program style to output in
        header: what to put in front of the cluster
        footer: what to put after the cluster
        ecp: capping ECP for the boundary region
        separate_pc: location for the separate point charge file (otherwise place in input file)
        """
        out = ''
        if 'header' in options:
            out += options['header']

        if options['program'] == 'orca':
            ecp = f' NewECP "{options["ecp"]}" end'

            form = '    {:<4}' + ' '*8 + ' {:> 13.8f}' * 3 + '\n'
            for atom, xyz, charge in self.qc_mol:
                out += form.format(atom, *xyz)

            form = '    {:<4}' + ' {:>7.4f}' + ' {:> 13.8f}' * 3
            for atom, xyz, charge in self.br_mol:
                out += form.format(atom + '>', charge, *xyz) + ecp + '\n'

            if 'separate_pc' in options:
                form = '{:>7.4f}' + ' {:> 13.8f}' * 3
                pc_out = f'{len(self.pc_mol)}\n'
                for atom, xyz, charge in self.pc_mol:
                    pc_out += form.format(charge, *xyz) + '\n'
                with open(options['separate_pc'], 'w') as f:
                    f.write(pc_out)
            else:
                for atom, xyz, charge in self.pc_mol:
                    out += form.format('Q', charge, *xyz) + '\n'

            out += '*'

        else:
            raise Exception(f'{options["program"]} is not yet supported')

        if 'footer' in options:
            out += options['footer']

        with open(outfile, 'w') as f:
            f.write(out)

    def recharged(self, charge):
        """ Change the net charge of non-qc region

        Distributes the increased or decreased charge between the atoms of the
        non-qc CMolecules

        :param charge: new charge
        """
        # Duplicate CMolecules
        qc_mol = CMolecule(self.qc_mol.geom)
        br_mol = CMolecule(self.br_mol.geom)
        pc_mol = CMolecule(self.pc_mol.geom)

        additional_charge = (charge - self.charge)/(len(br_mol) + len(pc_mol))

        br_mol.charges = br_mol.charges + additional_charge
        pc_mol.charges = pc_mol.charges + additional_charge

        return Cluster(qc_mol, br_mol, pc_mol)

    @staticmethod
    def from_radii(cmol, qc_radius, br_radius, center=(0, 0, 0)):
        """ Generates a cluster from regions defined by radii

        :param cmol: a CMolecule
        :params qc_radius, br_radius: radius of the given region
        :param center: center of the crystal (from which the radii radiate)
        """
        distances = np.linalg.norm(cmol.xyz - center, axis=1)

        cmols = []
        qc_region = distances <= qc_radius
        br_region = (qc_radius < distances) & (distances <= br_radius)
        pc_region = br_radius < distances
        for region in [qc_region, br_region, pc_region]:
            new_cmol = CMolecule()
            new_cmol.xyz = cmol.xyz[region]
            new_cmol.atoms = list(itertools.compress(cmol.atoms, region))
            new_cmol.charges = cmol.charges[region]
            cmols.append(new_cmol)

        return Cluster(*cmols)

    @staticmethod
    def from_rectangles(cmol, qc, br, center=(0, 0, 0)):
        """ Generates a cluster from regions defined by rectangular prisms

        :param cmol: a CMolecule
        :params qc, br: tuple(x, y, z) magnitude of cutoff (exclusive)
        :param center: center of the crystal from which the regions are generated

        +---+---+---+---+---+
        | * | * | * | * | * |
        +---+---+---+---+---+
        | * | b | b | b | * |
        +---+---+---+---+---+  * point charge region
        | * | b | q | b | * |  b boundary region
        +---+---+---+---+---+  q quantum region
        | * | b | b | b | * |
        +---+---+---+---+---+
        | * | * | * | * | * |
        +---+---+---+---+---+
        """
        center = np.array(center)

        def inside_indices(cmol, abc):

            def inside(xyz, abc):
                for coord, region in zip(xyz, abc):
                    if coord <= -region or region <= coord:
                        return False
                return True

            return {i for i, xyz in enumerate(cmol.xyz) if inside(xyz, abc)}

        qc_indices = inside_indices(cmol, qc - center)
        br_indices = inside_indices(cmol, br - center) - qc_indices
        pc_indices = set(range(len(cmol))) - qc_indices - br_indices

        cmols = []
        for indices in qc_indices, br_indices, pc_indices:
            new_cmol = CMolecule([data for i, data in enumerate(cmol) if i in indices])
            cmols.append(new_cmol)

        return Cluster(*cmols)
