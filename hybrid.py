import residue
import graph
import fep
import qscripts as q
import argparse
import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb

"""
Module that makes hybrid .pdb, .lib and .prm files of two given residues with a maximum common substructure
"""

parser = argparse.ArgumentParser(description="""
Command-line tool for creating a maximum common substructure hybrid residue 
""")
reqarg = parser.add_argument_group('Required')
reqarg.add_argument('res1', help='residue 1')
reqarg.add_argument('res2', help='residue 2')

args = parser.parse_args()


class Hybrid(object):
    def __init__(self, res1, res2):                 # reads in the .pdb, .lib, and .prm file for both residues
        self.res1 = residue.Residue(res1)
        self.res1.get_pdb(self.res1.pdb_file)
        self.res1.get_atoms(self.res1.lib)
        self.res1.get_bonds(self.res1.lib)
        self.res1.get_impropers(self.res1.lib)
        self.res1.get_connectivity(self.res1.mol)
        self.res1.get_atom_types(self.res1.ff_prm)
        self.res1.get_bond_types(self.res1.ff_prm)
        self.res1.get_angle_types(self.res1.ff_prm)
        self.res1.get_torsion_types(self.res1.ff_prm)
        self.res1.get_improper_types(self.res1.ff_prm)

        self.res2 = residue.Residue(res2)
        self.res2.get_pdb(self.res2.pdb_file)
        self.res2.get_atoms(self.res2.lib)
        self.res2.get_bonds(self.res2.lib)
        self.res2.get_impropers(self.res2.lib)
        self.res2.get_connectivity(self.res2.mol)
        self.res2.get_atom_types(self.res2.ff_prm)
        self.res2.get_bond_types(self.res2.ff_prm)
        self.res2.get_angle_types(self.res2.ff_prm)
        self.res2.get_torsion_types(self.res2.ff_prm)
        self.res2.get_improper_types(self.res2.ff_prm)

        self.protein_file = 'templates/protein/1bni_preproc_hbondopt_Q.pdb'
        self.mut_resn = 32
        self.mut_chaid_id = 'A'
        self.switching_atoms = {}
        self.ff_lib = 'templates/Qoplsaa.lib'

        self.hybrid_lib = 'tmp/hybrid.lib'
        self.hybrid_prm = 'tmp/hybrid.prm'
        self.hybrid_pdb = 'tmp/hybrid.pdb'
        self.hybrid_protein = 'tmp/protein_hybrid.pdb'


    def get_hybrid_atoms(self, res1, res2):
        """
        Function to find the atoms of the hybrid residue, returns these as a DataFrame
        """
        check = res1.atoms.atoms.isin(res2.atoms.atoms) # boolean check for atoms in res1.atoms are also in res2.atoms
        mcs = res1.atoms[check == True].copy()    # saves DataFrame with all atoms that are both in res1 and res2
        res1_rest = res1.atoms[check == False].copy()     # saves rest of atoms of res1 in DataFrame
        check = res2.atoms.atoms.isin(res1.atoms.atoms) # boolean check for atoms in res2.atoms are also in res1.atoms
        res2_rest = res2.atoms[check == False].copy()   # saves rest of atoms of res2 in DataFrame
        res2_rest.atoms = res2_rest.atoms.str.lower()   # makes res2 atoms lowercase
        self.hybrid_atoms = mcs._append([res1_rest, res2_rest], ignore_index=True)  # creates DataFrame with all hybrid atoms
        return self.hybrid_atoms

    def make_hybrid_pdb(self, pdb1, pdb2, hybrid_atoms):
        """
        Function to make a .pdb file of the hybrid residue
        """
        res1_check = pdb1.df['ATOM'].atom_name.isin(hybrid_atoms.atoms) # checks which atoms of pdb1 are in the hybrid_atoms DataFrane
        swap_atoms = hybrid_atoms.atoms.str.swapcase()   # Turns res2 atom_names to uppercase
        res2_check = pdb2.df['ATOM'].atom_name.isin(swap_atoms)     # checks which atoms of pdb2 are in the hybrid atoms DataFrame
        atoms_res1 = pdb1.df['ATOM'][res1_check == True].copy()     # takes hybrid atoms from pdb1
        atoms_res2 = pdb2.df['ATOM'][res2_check == True].copy()     # takes hybrid atoms from pdb2
        atoms_res2.atom_name = atoms_res2.atom_name.str.lower()     # makes res2 atoms lowercase
        hybrid_pdb = pd.concat([atoms_res1, atoms_res2], ignore_index=True)     # creates hybrid pdb DataFrame
        hybrid_pdb = hybrid_pdb.assign(atom_number=[x for x in range(1, len(hybrid_pdb)+1)],
                                       line_idx=[x for x in range(0, len(hybrid_pdb))],
                                       residue_name='HYB')    # renumbers hybrid pdb and
        hybrid_pdb.loc[:, 'chain_id'] = self.mut_chaid_id
        hybrid_pdb.loc[:, 'residue_number'] = self.mut_resn
        hybrid = PandasPdb()
        hybrid.df['ATOM'] = hybrid_pdb
        self.hybrid = hybrid.df['ATOM']
        hybrid.to_pdb(self.hybrid_pdb)
        return

    def make_hybrid_protein(self, protein):
        """
        Function to replace the WT residue with the hybrid residue in the protein.pdb file
        """
        self.protein = PandasPdb().read_pdb(protein)
        start = self.protein.df['ATOM'].loc[self.protein.df['ATOM']['residue_number'] < 32]
        end = self.protein.df['ATOM'].loc[self.protein.df['ATOM']['residue_number'] > 32]
        hybrid = pd.concat([start, self.hybrid, end], ignore_index=True)
        hybrid.loc[:, 'line_idx'] = [i for i in range(len(hybrid))]
        hybrid.loc[:, 'atom_number'] = [(i+1) for i in range(len(hybrid))]
        protein_hybrid = PandasPdb()
        protein_hybrid.df['ATOM'] = hybrid
        protein_hybrid.to_pdb(self.hybrid_protein)
        return

    def make_hybrid_lib(self):
        """
        Function to make a .lib file of the hybrid residue
        """

        # {HYB}
        with open('tmp/hybrid.lib', 'w') as lib_out:
            lib_out.write('{HYB}\n')

            # [atoms]
            lib_out.write('\n[atoms]\n')
            for index, row in self.hybrid_atoms.iterrows():
                lib_out.write(f"{index+1:>2} {row.atoms:>4}    {row.type:<4} {row.charge:>10}\n")         # prints the [atom] section of .lib

            # [bonds]
            lib_out.write('\n[bonds]\n')
            check = ~self.res2.bonds.isin(self.res1.bonds)  # inverted check which bonds are present in both residues
            self.res1.get_bond_types(self.res1.ff_prm)
            self.res2.get_bond_types(self.res2.ff_prm)
            self.hybrid_bonds = self.res1.bonds._append(self.res2.bonds[check.atom2 == True].copy(), ignore_index=True)
            self.hybrid_bonds.drop_duplicates(inplace=True, ignore_index=True)

            atoms1 = pd.concat([self.hybrid_bonds.atom1.isin(self.res1.atoms.atoms),
                               self.hybrid_bonds.atom1.isin(self.res2.atoms.atoms)],
                               axis=1, keys=['res1', 'res2'])
            atoms2 = pd.concat([self.hybrid_bonds.atom2.isin(self.res1.atoms.atoms),
                               self.hybrid_bonds.atom2.isin(self.res2.atoms.atoms)],
                               axis=1, keys=['res1', 'res2'])

            for i, bond in self.hybrid_bonds.iterrows():
                if atoms1.res1[i] == False and atoms1.res2[i] == True:
                    self.hybrid_bonds.at[i,'atom1'] = bond.atom1.lower()
                if atoms2.res1[i] == False and atoms2.res2[i] == True:
                    self.hybrid_bonds.at[i,'atom2'] = bond.atom2.lower()
                lib_out.write(f'{bond.atom1:>5} {bond.atom2:>5}\n')
                if (bond.atom1.isupper() and bond.atom2.islower()) or (bond.atom1.islower() and bond.atom2.isupper()):
                    if bond.atom1.isupper() and not bond.atom1 in self.switching_atoms:
                        self.switching_atoms[bond.atom1] = []
                    elif bond.atom2.isupper() and not bond.atom2 in self.switching_atoms:
                        self.switching_atoms[bond.atom2] = []
                    else:
                        pass
            for at in self.switching_atoms.keys():
                self.switching_atoms[at].append(self.res1.atoms.loc[self.res1.atoms.atoms == at].charge.values[0])
                self.switching_atoms[at].append(self.res2.atoms.loc[self.res2.atoms.atoms == at].charge.values[0])
                self.switching_atoms[at].append(self.res1.atoms.loc[self.res1.atoms.atoms == at].type.values[0])
                self.switching_atoms[at].append(self.res2.atoms.loc[self.res2.atoms.atoms == at].type.values[0])

            # [connections]
            lib_out.write('\n[connections]\n')
            lib_out.write(' head N\n')
            lib_out.write(' tail C\n')

            # [impropers]
            lib_out.write('\n[impropers]\n')
            if self.res1.get_impropers(self.res1.lib).empty and self.res2.get_impropers(self.res2.lib).empty:   # skip if no impropers in both residues
                pass
            else:
                self.res1.get_impropers(self.res1.lib)  # get the impropers from res1
                if not self.res1.impropers.empty:   # if there are
                    for i, row in self.res1.impropers.iterrows():   # add them
                        lib_out.write(f'{row.atom1:>5} {row.atom2:>5} {row.atom3:>5} {row.atom4:>5}\n')
                self.res2.get_impropers(self.res2.lib)  # get the impropers from res2
                if not self.res2.impropers.empty:   # if there are
                    check = self.res2.impropers.isin(self.res1.impropers)
                    for i, row in self.res2.impropers.iterrows():   # add them
                        if check.loc[i,:].all() == True:
                            continue
                        else:
                            if check.atom1[i] == True:
                                lib_out.write(f'{row.atom1:>5} ')
                            else:
                                lib_out.write(f'{row.atom1.lower():>5} ')
                            if check.atom2[i] == True:
                                lib_out.write(f'{row.atom2:>5} ')
                            else:
                                lib_out.write(f'{row.atom2.lower():>5} ')
                            if check.atom3[i] == True:
                                lib_out.write(f'{row.atom3:>5} ')
                            else:
                                lib_out.write(f'{row.atom3.lower():>5} ')
                            if check.atom4[i] == True:
                                lib_out.write(f'{row.atom4:>5}\n')
                            else:
                                lib_out.write(f'{row.atom5.lower():>5}\n')

            # [charge_groups]
            lib_out.write('\n[charge_groups]\n')

            lib_out.close()
            return

    def atom2type(self, atom):
        """
        Function that takes an atom name and returns an atom type
        """
        type = self.hybrid_atoms[self.hybrid_atoms.atoms == atom].type.values[0]
        return type

    def add_hybrid_prms(self):
        """
        Function to add parameters of the hybrid residue to the force field .prm file
        """

        hybrid_prms = {'atoms': [], 'bonds': [], 'angles': [], 'torsions': [], 'impropers': []}
        hybrid_type = {'type': None, 'A1': '0.00', 'A2': '0.00', 'B1': '0.00', 'A3': '0.00', 'B2': '0.00', 'mass': '1.008'}
        missing_atoms = list(set(self.hybrid_atoms.type.values).difference(self.res1.ff_atoms.type.values))
        if missing_atoms:
            for at_type in missing_atoms:
                hybrid_type['type'] = at_type
                hybrid_prms['atoms'].append(hybrid_type)
                # self.res1.ff_atoms = self.res1.ff_atoms._append(hybrid_type, ignore_index=True)

        hybrid_bond = {'type1': None, 'type2': None, 'kb': '0.0', 'r0': '1.090'}
        bonds = self.hybrid_bonds.applymap(self.atom2type)
        for b_type in bonds[['atom1', 'atom2']].values.tolist():
            if b_type in self.res1.ff_bonds[['type1', 'type2']].values.tolist():
                continue
            elif b_type[::-1] in self.res1.ff_bonds[['type1', 'type2']].values.tolist():
                continue
            else:
                hybrid_bond['type1'] = b_type[0]
                hybrid_bond['type2'] = b_type[1]
                hybrid_prms['bonds'].append(hybrid_bond.copy())
                # self.res1.ff_bonds = self.res1.ff_bonds._append(hybrid_bond, ignore_index=True)

        hybrid_angle = {'type1': None, 'type2': None, 'type3': None, 'kth': '0.0', 'th0': '0.000'}
        for a_type in self.hybrid_angles:
            if a_type in self.res1.ff_angles[['type1', 'type2', 'type3']].values.tolist():
                continue
            elif a_type[::-1] in self.res1.ff_angles[['type1', 'type2', 'type3']].values.tolist():
                continue
            else:
                hybrid_angle['type1'] = a_type[0]
                hybrid_angle['type2'] = a_type[1]
                hybrid_angle['type3'] = a_type[2]
                hybrid_prms['angles'].append(hybrid_angle.copy())
                # self.res1.ff_angles = self.res1.ff_angles._append(hybrid_angle, ignore_index=True)

        hybrid_torsion = {'type1': None, 'type2': None, 'type3': None, 'type4': None,
                           'kph': '0.000', 'n': '1.000', 'd': '0.000', 'p': '1.00'}
        for t_type in self.hybrid_torsions:
            if t_type in self.res1.ff_torsions[['type1', 'type2', 'type3', 'type4']].values.tolist():
                continue
            elif t_type[::-1] in self.res1.ff_torsions[['type1', 'type2', 'type3', 'type4']].values.tolist():
                continue
            else:
                hybrid_torsion['type1'] = t_type[0]
                hybrid_torsion['type2'] = t_type[1]
                hybrid_torsion['type3'] = t_type[2]
                hybrid_torsion['type4'] = t_type[3]
                hybrid_prms['torsions'].append(hybrid_torsion.copy())
                # self.res1.ff_torsions = self.res1.ff_torsions._append(hybrid_torsion, ignore_index=True)

        with open('templates/Qoplsaa.prm', 'r') as prm_in, open('tmp/hybrid.prm', 'w') as prm_out:
            for line in prm_in:
                if line.startswith('!Hybrid'):
                    prm_out.write(line)
                    for prm in hybrid_prms[line.split('-')[2]]:
                        prm_out.write(('{:<8} '*len(prm.values())+'\n').format(*prm.values()))
                    continue
                prm_out.write(line)

    def make_graph(self):
        graph.Graph(self)
        return

    def qprep(self):
        cyx = ['at1', 'at2']

        self.center = '33.963 31.279 29.674'
        self.radius = '25'
        self.solvent = 'HOH'
        q.qprep(self.hybrid_lib, self.ff_lib, self.hybrid_prm, self.hybrid_protein,
                self.center, self.radius, self.solvent, cyx)
        return

    def make_fep_file(self):
        fep.Fep(self).execute()


test = Hybrid(res1=args.res1, res2=args.res2)
test.get_hybrid_atoms(test.res1, test.res2)
test.make_hybrid_pdb(test.res1.pdb, test.res2.pdb, test.hybrid_atoms)
test.make_hybrid_lib()
test.make_graph()
test.add_hybrid_prms()
test.make_hybrid_protein(test.protein_file)
test.qprep()
test.make_fep_file()

""" softcode:
make sure that the Qoplsaa.lib file gets read for the residue library
make sure that the Qoplsaa.prm file gets read for the residue parameters
have a template pdb for each amino acid"""