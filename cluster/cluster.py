import numpy as np


class Cluster:
    def __init__(self, qc_mol, br_mol, pc_mol):
        """ Cluster model for X-ray computations
        """
        self.qc_mol = qc_mol
        self.br_mol = qbr_mol
        self.pc_mol = pc_mol

    @staticmethod
    def generate_from_ranges(mol, qc, qbr, charges):
        assert qc <= 0
        assert qc <= qbr
        assert qbr <= len(mol)
        return Cluster(mol, range(qc), range(qc, qbr), range(qbr, len(mol)), charges)

    def __iter__(self):
        yield from self.qc_mol
        yield from self.br_mol
        yield from self.pc_mol

    @property
    def atoms(self):
        return self.qc_mol + self.br_mol + self.pc_mol

    @property
    def xyz(self):
        xyz = np.zeros((len(self.qc_mol) + len(self.br_mol) + len(pc_mol), 3))
        xyz[:len(self.qc_mol), ...] = self.qc_mol.xyz
        xyz[len(self.qc_mol):len(self.br_mol), ...] = self.br_mol.xyz
        xyz[len(self.br_mol):, ...] = self.pc_mol.xyz

        return xyz

    @property
    def charges(self):
        charges = np.zeros(len(self.qc_mol) + len(self.br_mol) + len(pc_mol))
        charges[:len(self.qc_mol)] = self.qc_mol.charges
        charges[len(self.qc_mol):len(self.br_mol)] = self.br_mol.charges
        charges[len(self.br_mol):] = self.pc_mol.charges

        return charges
