import os
import sys
import pandas as pd
from biopandas.pdb import PandasPdb
from openbabel import openbabel as ob
from openbabel import pybel

"""
Module for conversion of molecule .pdb, .lib and .prm files to data structures.
"""


class Residue(object):
    def __init__(self, name):
        self.name = name
        self.pdb_file = f'templates/residues/{name}/{name}.pdb'
        self.lib = f'templates/residues/{name}/{name}.lib'
        # self.prm = f'templates/residues/{name}/{name}.prm'
        self.mol = next(pybel.readfile('pdb', self.pdb_file))

    def read_file(self, file, start_cond, end_cond, *args):
        """
        Function for reading .lib and .prm files, returns list of lines
        """
        if args:
            end = [end_cond] + [arg for arg in args]
        else:
            end = [end_cond]
        read = False
        lines = []
        with open(file, 'r') as read_file:
            for line in read_file:
                if any(cond in line for cond in end):
                    read = False
                if read:
                    if not line.strip():
                        continue
                    else:
                        lines.append(line.split('#')[0].split())
                if start_cond in line:
                    read = True
        return lines

    def get_pdb(self, pdb_file):
        """
        Function for reading a .pdb file with Biopandas, returns a Biopandas DataFrame
        """
        self.pdb = PandasPdb().read_pdb(pdb_file)
        return self.pdb

    def atom2type(self, atom):
        """
        Function that takes an atom name and returns an atom type
        """
        type = self.atoms[self.atoms.atoms == atom].type.values[0]
        return type

    def type2atom(self, type):
        """
        Function that takes an atom type and returns an atom name
        """
        atom = self.atoms[self.atoms.type == type].atoms.values[0]
        return atom

    def get_atoms(self, lib):
        """
        Function that returns [atom] entries from a .lib file in a DataFrame object
        """
        atoms = self.read_file(lib, '[atoms]', '[bonds]')
        self.atoms = pd.DataFrame(atoms, columns=['index', 'atoms', 'type', 'charge']).set_index('index')
        return self.atoms

    def get_atom_types(self, prm):
        """
        Function that returns [atom_types] entries from a .prm file in a DataFrame object
        """
        atom_types = self.read_file(prm, '[atom_types]', '[bonds]')
        atom_types = pd.DataFrame(atom_types, columns=['type', 'A1', 'A2', 'B1', 'A3', 'B2', 'mass'])
        self.atoms = pd.merge(self.atoms, atom_types, on='type')
        return self.atoms

    def get_bonds(self, lib):
        """
        Function that returns [bond] entries from a .lib file in a DataFrame object
        """
        bonds = self.read_file(lib, '[bonds]', '[connections]', '[charge_groups]', '[impropers]')
        self.bonds = pd.DataFrame(bonds, columns=['atom1', 'atom2'])
        return self.bonds

    # def get_bond_types(self, prm):
    #     """
    #     Function that returns [bond_types] entries from a .prm file in a DataFrame object
    #     """
    #     bond_types = self.read_file(prm, '[bonds]', '[angles]')
    #     bond_types = pd.DataFrame(bond_types, columns=['type1', 'type2', 'kb', 'r0'])
    #     bond_types.type1 = bond_types.type1.map(self.type2atom)
    #     bond_types.type2 = bond_types.type2.map(self.type2atom)
    #     self.bonds = pd.merge(self.bonds, bond_types, how='right', left_on=['atom1', 'atom2'], right_on=['type1', 'type2'])
    #     self.bonds.atom1 = self.bonds.type1
    #     self.bonds.atom2 = self.bonds.type2
    #     self.bonds.type1 = self.bonds.type1.map(self.atom2type)
    #     self.bonds.type2 = self.bonds.type2.map(self.atom2type)
    #     # self.bonds.insert(loc=0, column='type1', value=bond_types['type1'].map(self.atom2type))
    #     # self.bonds.insert(loc=2, column='type2', value=bond_types['type2'].map(self.atom2type))
    #     return self.bonds

    # def get_angles(self, prm):
    #     """
    #     Function that returns [angles] and [angle_types] entries from a .prm file in a DataFrame object
    #     """
    #     angle_types = self.read_file(prm, '[angles]', '[torsions]')
    #     angle_types = pd.DataFrame(angle_types, columns=['type1', 'type2', 'type3', 'kth', 'th0'])
    #     angle_types.insert(loc=0, column='atom1', value=angle_types['type1'].map(self.type2atom))
    #     angle_types.insert(loc=2, column='atom2', value=angle_types['type2'].map(self.type2atom))
    #     angle_types.insert(loc=4, column='atom3', value=angle_types['type3'].map(self.type2atom))
    #     self.angles = angle_types
    #     return self.angles

    # def get_torsions(self, prm):
    #     """
    #     Function that returns [torsions] and [torsion_types] entries from a .prm file in a DataFrame object
    #     """
    #     torsion_types = self.read_file(prm, '[torsions]', '[impropers]')
    #     torsion_types = pd.DataFrame(torsion_types, columns=['type1', 'type2', 'type3', 'type4', 'Kph', 'n', 'd', 'p'])
    #     torsion_types.insert(loc=0, column='atom1', value=torsion_types['type1'].map(self.type2atom))
    #     torsion_types.insert(loc=2, column='atom2', value=torsion_types['type2'].map(self.type2atom))
    #     torsion_types.insert(loc=4, column='atom3', value=torsion_types['type3'].map(self.type2atom))
    #     torsion_types.insert(loc=6, column='atom4', value=torsion_types['type4'].map(self.type2atom))
    #     self.torsions = torsion_types
    #     return self.torsions

    # def get_impropers(self, lib):
    #     """
    #     Function that returns [impropers] entries from a .lib file in a DataFrame object
    #     """
    #     impropers = self.read_file(lib, '[impropers]', '[charge_groups]', 'EOF')
    #     impropers = pd.DataFrame(impropers, columns=['atom1', 'atom2', 'atom3', 'atom4'])
    #     self.impropers = impropers
    #     print(self.impropers)
    #     return self.impropers

    # def get_impropers(self, prm):
    #     """
    #     Function that returns [impropers] and [improper_types] entries from a .prm file in a DataFrame object
    #     """
    #     improper_types = self.read_file(prm, '[impropers]', 'EOF')
    #     improper_types = pd.DataFrame(improper_types, columns=['type1', 'type2', 'type3', 'type4', 'kx', 'x0'])
    #     improper_types.insert(loc=0, column='atom1', value=improper_types['type1'].map(self.type2atom))
    #     improper_types.insert(loc=2, column='atom2', value=improper_types['type2'].map(self.type2atom))
    #     improper_types.insert(loc=4, column='atom3', value=improper_types['type3'].map(self.type2atom))
    #     improper_types.insert(loc=6, column='atom4', value=improper_types['type4'].map(self.type2atom))
    #     self.impropers = improper_types
    #     return self.impropers

    def get_charge_groups(self, lib):
        """
        Function that returns [charge_groups] entries from a .lib file in a DataFrame object
        """
        charge_groups = self.read_file(lib, '[charge_groups', '*-', 'EOF')
        self.charge_groups = pd.DataFrame(charge_groups, columns=[f'atom{i}' for i in range(max(len(j) for j in charge_groups))])
        return self.charge_groups

    def get_connectivity(self, mol):
        """
        Function that returns all angles in a given residue
        """

        self.angles = []
        for angle in ob.OBMolAngleIter(mol.OBMol):
            atoms = [self.atoms.atoms[atom] for atom in angle]
            atoms[0], atoms[1] = atoms[1], atoms[0]
            self.angles.append(atoms)

        self.angle_types = []
        for angle in self.angles:
            types = [self.atom2type(atom) for atom in angle]
            self.angle_types.append(types)

        self.torsions = []
        for torsion in ob.OBMolTorsionIter(mol.OBMol):
            atoms = [self.atoms.atoms[atom] for atom in torsion]
            self.torsions.append(atoms)
        self.torsion_types = []
        for torsion in self.torsions:
            types = [self.atom2type(atom) for atom in torsion]
            self.torsion_types.append(types)

        return self.angles, self.torsions, self.angle_types, self.torsion_types


test = Residue('ALA')
test.get_pdb(test.pdb_file)
test.get_atoms(test.lib)
print(test.atoms)
# test.get_atom_types(test.prm)
test.get_bonds(test.lib)
print(test.bonds)
test.get_connectivity(test.mol)
# print(test.__dict__)
# test.get_bond_types(test.prm)

# test2 = Residue('ARG/n-propylguanidine')
# test2.get_atoms(test2.lib)
# # test2.get_atom_types(test2.prm)
# # test2.get_bonds(test2.lib)
# # test2.get_bond_types(test2.prm)
# # test2.get_angles(test2.prm)
# # test2.get_torsions(test2.prm)
# test2.get_impropers(test2.prm)
# test2.get_charge_groups(test2.lib)
