from difflib import restore
from operator import truediv
from os import truncate
from xml.etree.ElementTree import TreeBuilder
import residue
import argparse
import pandas as pd 
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
        self.res1 = residue.Residue(res1)           # stores all this data in the object self.res
        self.res2 = residue.Residue(res2)           # self.res.atoms is DataFrame with all values
        self.res1.get_atoms(self.res1.lib)
        self.res1.get_atom_types(self.res1.prm)
        self.res2.get_atoms(self.res2.lib)
        self.res2.get_atom_types(self.res2.prm)
        self.res1.get_pdb(self.res1.pdb_file)
        self.res2.get_pdb(self.res2.pdb_file)
        self.res1.get_bonds(self.res1.lib)
        self.res2.get_bonds(self.res2.lib)
        self.res1.get_bond_types1(self.res1.lib)
        self.res2.get_bond_types1(self.res2.lib)
        self.res1.get_bond_types2(self.res1.prm)
        self.res2.get_bond_types2(self.res2.prm)
        self.res1.get_angles(self.res1.prm)
        self.res2.get_angles(self.res2.prm)
        self.res1.get_torsions(self.res1.prm)
        self.res2.get_torsions(self.res2.prm)
        #print(self.res1.atoms)


    def get_hybrid_atoms(self):
        """
        Function to find the atoms of the hybrid residue, returns these as a DataFrame
        """
        check_res1 = self.res1.atoms['atoms'].isin(self.res2.atoms['atoms'])
        check_res2 = self.res2.atoms['atoms'].isin(self.res1.atoms['atoms'])
        mcs = self.res1.atoms[check_res1 == True].copy()    # saves DataFrame with all atoms that are both in res1 and res2
        res1_rest = self.res1.atoms[check_res1 == False].copy()     # saves rest of atoms of res1 in DataFrame
        res2_rest = self.res2.atoms[check_res2 == False].copy()     # saves rest of atoms of res2 in Dataframe
        res2_rest.loc[:, 'atoms'] = res2_rest.loc[:, 'atoms'].str.lower()  # makes res2 atoms lowercase
        self.hybrid_atoms = mcs.append([res1_rest, res2_rest], ignore_index=True)   # creates DataFrame with all hybrid atoms
        #print(self.hybrid_atoms)
        return self.hybrid_atoms

    def make_hybrid_pdb(self, pdb1, pdb2, hybrid_atoms):
        """
        Function to make a .pdb file of the hybrid residue
        """
        pdb1_atoms = pdb1.df['ATOM']
        res1_check = pdb1_atoms['atom_name'].isin(hybrid_atoms['atoms'])    # checks which atoms of pdb1 are in the hybrid_atoms DataFrame  
        swap_atoms = hybrid_atoms['atoms'].str.swapcase()   # Turns res2 atom_names to uppercase and res1 to lowercase
        compare_atoms = pdb2.df['ATOM'].atom_name       # makes DataFrame of atom_names of pdb2
        res2_check = compare_atoms.isin(swap_atoms)     # checks which atoms of pdb2 are in the hybrid atoms DataFrame
        atoms_res1 = pdb1.df['ATOM'][res1_check == True].copy()     # takes hybrid atoms from pdb1
        atoms_res2 = pdb2.df['ATOM'][res2_check == True].copy()     # takes hybrid atoms from pdb2
        atoms_res2.loc[:, 'atom_name'] = atoms_res2.loc[:, 'atom_name'].str.lower()     # makes res2 atoms lowercase
        hybrid_pdb = atoms_res1.append(atoms_res2, ignore_index=True)       # creates hybrid pdb DataFrame
        hybrid_pdb = hybrid_pdb.assign(atom_number=[x for x in range(1, len(hybrid_pdb)+1)], residue_name='HYB')    # renumbers hybrid pdb and
        hybrid = PandasPdb()
        hybrid.df['ATOM'] = hybrid_pdb
        hybrid.to_pdb('tmp/hybrid.pdb')
        #print(hybrid_pdb)
        return

    def make_hybrid_lib(self, lib1, lib2):
        """
        Function to make a .lib file of the hybrid residue
        """

        # {HYB}
        with open('tmp/hybrid.lib', 'w') as lib_out:
            lib_out.write('{HYB}\n')

            # [atoms]
            lib_out.write('\n[atoms]\n')
            for index, row in self.hybrid_atoms.iterrows():
                lib_out.write(f"{index+1:>2} {row.atoms:>4}    hyb.{row.atoms:<4} {row.charge:>10}\n")         # prints the [atom] section of .lib
            # [bonds]
            lib_out.write('\n[bonds]\n')
            merged_data = pd.merge(self.res1.bonds, self.res2.bonds, on=['atom1', 'atom2'], how='outer', indicator=True) #check which bonds are present in both residues
            merged_data['_merge'] = merged_data['_merge'].replace({'both': True, 'right_only': False,'left_only': True}) #True of false for left_only??
            merged_data.loc[merged_data['_merge'] == False, ['atom1', 'atom2']] = merged_data.loc[merged_data['_merge'] == False, ['atom1', 'atom2']].apply(lambda x: x.str.lower())
            self.hybrid_bonds = merged_data.drop('_merge', axis=1)
            for i, row in self.hybrid_bonds.iterrows():
                lib_out.write(f'{row.atom1:>5} {row.atom2:>5}\n')
            # [impropers]
            lib_out.write('\n[impropers]\n')
            merged_data = pd.merge(self.res1.get_impropers(self.res1.lib), self.res2.get_impropers(self.res2.lib), on=['atom1', 'atom2', 'atom3', 'atom4'], how='outer', indicator=True) #check which impropers are present in both residues
            merged_data['_merge'] = merged_data['_merge'].replace({'both': True, 'right_only': False,'left_only': False})
            merged_data.loc[merged_data['_merge'] == False, ['atom1', 'atom2','atom3', 'atom4']] = merged_data.loc[merged_data['_merge'] == False, ['atom1', 'atom2', 'atom3', 'atom4']].apply(lambda x: x.str.lower())
            self.hybrid_impropers = merged_data.drop('_merge', axis=1)
            for i, row in self.hybrid_impropers.iterrows():
                lib_out.write(f'{row.atom1:>5} {row.atom2:>5} {row.atom3:>5} {row.atom4:>5}\n')
            # [charge_groups]
            lib_out.write('\n[charge_groups]\n')
            column_names_res1 = self.res1.get_charge_groups(self.res1.lib).columns.tolist()
            column_names_res2 = self.res2.get_charge_groups(self.res2.lib).columns.tolist()
            merged_data = pd.merge(self.res1.get_charge_groups(self.res1.lib), self.res2.get_charge_groups(self.res2.lib), on=column_names_res1, how='outer', indicator=True)
            merged_data['_merge'] = merged_data['_merge'].replace({'both': True, 'right_only': False,'left_only': True}) 
            for index, row in merged_data.loc[merged_data['_merge'] == False].iterrows():
                atom_columns = [column for column in row.index if column.startswith('atom')]
                merged_data.loc[index, atom_columns] = row[atom_columns].str.lower()
            self.hybrid_impropers = merged_data.drop('_merge', axis=1)
            for i, row in self.hybrid_impropers.iterrows():
                group_atoms = ' '.join(f'{str(atom):<5}' if pd.notnull(atom) else '     ' for atom in row.values)
                lib_out.write(f'{group_atoms}\n')
            lib_out.close()
            return

    def make_hybrid_prm(self, prm1, prm2):
        """
        Function to make a .prm file of the hybrid residue
        """
        with open('tmp/hybrid.prm', 'w') as prm_out:
            prm_out.write('{HYB}\n')
            prm_out.write('[options]\n')
            # [atoms]
            prm_out.write('\n[atom_types]\n')
            self.res1.atom_types.rename(columns={self.res1.atom_types.columns[5]: 'Bvdw2'}, inplace=True) #Bvdw2 is actually Bvdw2&3 but it has difficulties reading the & character. 
            self.res2.atom_types.rename(columns={self.res2.atom_types.columns[5]: 'Bvdw2'}, inplace=True)
            merged_data = pd.merge(self.res1.atom_types, self.res2.atom_types, on= ['type', 'Avdw1', 'Avdw2', 'Bvdw1', 'Avdw3', 'Bvdw2','mass','atoms'], how='outer', indicator=True) #check which atom types are present in both residues
            merged_data['_merge'] = merged_data['_merge'].replace({'both': True, 'right_only': False,'left_only': True})
            merged_data.loc[merged_data['_merge'] == False, ['type', 'Avdw1', 'Avdw2', 'Bvdw1', 'Avdw3', 'Bvdw2','mass','atoms']] = merged_data.loc[merged_data['_merge'] == False, ['type', 'Avdw1', 'Avdw2', 'Bvdw1', 'Avdw3', 'Bvdw2','mass']].apply(lambda x: x.str.lower())
            self.hybrid_atom_types = merged_data.drop(['_merge','atoms'], axis=1)
            for i, row in self.hybrid_atom_types.iterrows():
                prm_out.write(f"{row.type:>7} {row.Avdw1:>7} {row.Avdw2:>7} {row.Bvdw1:>7} {row.Avdw3:>7} {row.Bvdw2:>7} {row.mass:>7}\n")
            #[bonds]
            prm_out.write('\n[bonds]\n')
            #Identify hybrid bond types, bond types only present in residue1, and residues only present in residue2
            self.res1.bond_types['merge'] = 'left_only'
            self.res2.bond_types['merge'] = 'right_only'
            for i i