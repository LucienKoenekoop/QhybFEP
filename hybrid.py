import residue
import argparse
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
        # self.res1.get_bond_types(self.res1.prm)
        # self.res2.get_bond_types(self.res2.prm)
        self.res2.get_bonds(self.res2.lib)
        self.res1.get_angles(self.res1.prm)
        self.res2.get_angles(self.res2.prm)
        self.res1.get_torsions(self.res1.prm)
        self.res2.get_torsions(self.res2.prm)


    def get_hybrid_atoms(self):
        """
        Function to find the atoms of the hybrid residue, returns these as a DataFrame
        """
        check = self.res1.atoms.isin({'atoms': self.res2.atoms['atoms']})['atoms']  # boolean check for atoms in self.res1.atoms are also in self.res2.atoms
        mcs = self.res1.atoms[check == True].copy()    # saves DataFrame with all atoms that are both in res1 and res2
        res1_rest = self.res1.atoms[check == False].copy()     # saves rest of atoms of res1 in DataFrame
        check, other = check.align(self.res2.atoms, axis=0)
        res2_rest = self.res2.atoms[check.loc[0:len(self.res2.atoms)-1] == False].copy()     # saves rest of atoms of res2 in DataFrame
        res2_rest.loc[:, 'atoms'] = res2_rest.loc[:, 'atoms'].str.lower()  # makes res2 atoms lowercase
        self.hybrid_atoms = mcs.append([res1_rest, res2_rest], ignore_index=True)       # creates DataFrame with all hybrid atoms
        return self.hybrid_atoms

    def make_hybrid_pdb(self, pdb1, pdb2, hybrid_atoms):
        """
        Function to make a .pdb file of the hybrid residue
        """
        res1_check = pdb1.df['ATOM'].isin({'atom_name': hybrid_atoms['atoms']})['atom_name']    # checks which atoms of pdb1 are in the hybrid_atoms DataFrane
        swap_atoms = hybrid_atoms['atoms'].str.swapcase()   # Turns res2 atom_names to uppercase
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
            check = ~self.res2.bonds.isin(self.res1.bonds)  # inverted check which bonds are present in both residues
            print(check)
            self.res1.get_bond_types(self.res1.prm)
            self.res2.get_bond_types(self.res2.prm)
            print(self.res1.bonds)
            print(self.res2.bonds)
            # self.hybrid_bonds = self.res1.bonds.append(self.res2.bonds[check.atom2.loc[0:len(self.res2.bonds)-1] == True])  # creates hybrid bonds DataFrame
            self.hybrid_bonds = self.res1.bonds.append(self.res2.bonds[check.atom2 == True].copy())
            print(self.hybrid_bonds)
            a = 0
            for i, row in self.hybrid_bonds.iterrows():
                if int(i) == a:
                    lib_out.write(f'{row.atom1:>5} {row.atom2:>5}\n')
                else:
                    if int(i) == (a-1):
                        lib_out.write(f'{row.atom1.lower():>5} {row.atom2.lower():>5}\n')
                    else:
                        lib_out.write(f'{row.atom1:>5} {row.atom2.lower():>5}\n')
                        a = int(i)+1
                a += 1

            # [connections]
            # lib_out.write('\n[connections]\n')
            # lib_out.write(' head N')
            # lib_out.write(' tail C')

            # [impropers]
            lib_out.write('\n[impropers]\n')
            if self.res1.get_impropers(self.res1.prm).empty and self.res2.get_impropers(self.res2.prm).empty:   # skip if no impropers in both residues
                pass
            else:
                self.res1.get_impropers(self.res1.prm)  # get the impropers from res1
                if not self.res1.impropers.empty:   # if there are
                    for i, row in self.res1.impropers.iterrows():   # add them
                        lib_out.write(f'{row.atom1:>5} {row.atom2:>5} {row.atom3:>5} {row.atom4:>5}\n')
                self.res2.get_impropers(self.res2.prm)  # get the impropers from res2
                if not self.res2.impropers.empty:   # if there are
                    for i, row in self.res2.impropers.iterrows():   # add them
                        lib_out.write(f'{row.atom1.lower():>5} {row.atom2.lower():>5} {row.atom3.lower():>5} {row.atom4.lower():>5}\n')

            # [charge_groups]
            lib_out.write('\n[charge_groups]\n')
            # self.res1.get_charge_groups(lib1)
            # self.res2.get_charge_groups(lib2)
            # for i, item in self.res2.charge_groups.iteritems():
            #     if item.values[0] in self.res1.charge_groups.values:
            #         pass
            #         # print(item.values[0])
            #         # print('True')
            #     else:
            #         item.values[0] = item.values[0].lower()
            # print(self.res1.charge_groups)
            # print(self.res2.charge_groups)

            lib_out.close()
            return

    def make_hybrid_prm(self, prm1, prm2):
        """
        Function to make a .prm file of the hybrid residue
        """

        with open('tmp/hybrid.prm', 'w') as prm_out:
            prm_out.write('[options]\n')
            prm_out.write('\n[atom_types]\n')
            for i, atom in self.res1.atoms.iterrows():
                prm_out.write(f'hyb.{atom.atoms:<8} {atom.A1:>8} {atom.A2:>10} {atom.B1:>10} {atom.A3:>10} {atom.B2:>10} {atom.mass:>10}\n')
            for i, atom in self.res2.atoms.iterrows():
                if atom.atoms in self.res1.atoms.atoms.values:
                    pass
                else:
                    prm_out.write(f'hyb.{atom.atoms.lower():<8} {atom.A1:>8} {atom.A2:>10} {atom.B1:>10} {atom.A3:>10} {atom.B2:>10} {atom.mass:>10}\n')

            prm_out.write('\n[bonds]\n')
            a = 0
            # print(self.hybrid_bonds)
            for i, row in self.hybrid_bonds.iterrows():
                if int(i) == a:
                    prm_out.write(f'hyb.{row.atom1:<5} hyb.{row.atom2:<5}\n')
                else:
                    if int(i) == (a-1):
                        prm_out.write(f'hyb.{row.atom1.lower():<5} hyb.{row.atom2.lower():<5}\n')
                    else:
                        prm_out.write(f'hyb.{row.atom1:<5} hyb.{row.atom2.lower():<5}\n')
                        a = int(i)+1
                a += 1

            # print(self.res1.bonds)
            # print(self.res1.angles)
            # print(self.res1.torsions)
            # print(self.res1.impropers)


test = Hybrid(res1=args.res1, res2=args.res2)
test.get_hybrid_atoms()
test.make_hybrid_pdb(test.res1.pdb, test.res2.pdb, test.hybrid_atoms)
test.make_hybrid_lib(test.res1.lib, test.res2.lib)
test.make_hybrid_prm(test.res1.prm, test.res2.prm)

""" softcode:
make sure that the Qoplsaa.lib file gets read for the residue library
make sure that the Qoplsaa.prm file gets read for the residue parameters
have a template pdb for each amino acid"""