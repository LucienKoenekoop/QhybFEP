import pandas as pd
class Fep(object):

    def __init__(self, hybrid):
        self.hybrid = hybrid
        self.res1 = hybrid.res1
        self.res2 = hybrid.res2
        self.backbone = ['N', 'H', 'CA', 'HA', 'C', 'O']
        self.switching_atoms = hybrid.switching_atoms
        # self.atom_types = [['DUM', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '1.0080']]
        self.mutation = 'A2V'
        return

    def get_top(self):
        q_atoms = []
        with open('tmp/protein_hybrid.pdb', 'r') as topology:
            for atom in topology:
                if 'HYB' in atom:
                    q_atoms.append((atom.split()[1], atom.split()[2]))
            topology.close()
        return q_atoms

    def charges(self):
        charges = []
        for i, charge in self.hybrid.hybrid_atoms.iterrows():
            if charge.atoms in self.backbone:
                charges.append((charge.charge, charge.charge))
            elif charge.atoms in self.switching_atoms:
                charges.append((self.switching_atoms[charge.atoms][0], self.switching_atoms[charge.atoms][1]))
            elif charge.atoms.isupper():
                charges.append((charge.charge, '0.0000'))
            elif charge.atoms.islower():
                charges.append(('0.0000', charge.charge))
        return charges

    def atom_types(self):
        atom_types = [['DUM', '0.00', '0.00', '0.00', '0.00', '0.00', '1.008']]
        for i, type in self.hybrid.hybrid_atoms.iterrows():
            if not type.type in [item[0] for item in atom_types]:
                atom_types.append(self.res1.ff_atoms.loc[self.res1.ff_atoms['type'] == type.type].values.tolist()[0])
        return atom_types

    def change_atoms(self):
        types = []
        for i, type in self.hybrid.hybrid_atoms.iterrows():
            if type.atoms in self.backbone:
                types.append((type.type, type.type))
            elif type.atoms in self.switching_atoms:
                types.append((self.switching_atoms[type.atoms][2], self.switching_atoms[type.atoms][3]))
            elif type.atoms.isupper():
                types.append((type.type, 'DUM'))
            elif type.atoms.islower():
                types.append(('DUM', type.type))
        return types

    def excluded_pairs(self, q_atoms):
        pairs = []
        res1 = self.hybrid.hybrid_atoms.atoms.loc[self.hybrid.hybrid_atoms['atoms'].str.isupper()].values.tolist()
        res2 = self.hybrid.hybrid_atoms.atoms.loc[self.hybrid.hybrid_atoms['atoms'].str.islower()].values.tolist()
        for atom in self.backbone:
            if atom in res1:
                res1.remove(atom)
        for at_i in res1:
            for at_j in res2:
                pairs.append((self.q_atnum(at_i, q_atoms), self.q_atnum(at_j, q_atoms)))
        return pairs

    def softcore(self):
        softcore = []
        for i, atom in self.hybrid.hybrid_atoms.iterrows():
            if atom.atoms in self.backbone or atom.atoms in self.switching_atoms:
                softcore.append(('0', '0'))
            elif atom.atoms.isupper():
                softcore.append(('0', '20'))
            elif atom.atoms.islower():
                softcore.append(('20', '0'))
        return softcore

    # def bond_types(self):
    def angle_types(self):
        angle_types = []
        return
    def zero_angles(self, q_atoms):
        zero_angles = []
        for angle in self.hybrid.fep_angles:
            forward = [x.upper() for x in angle]
            reverse = forward[::-1]
            if not forward in self.res1.angles.values.tolist() \
            and not reverse in self.res1.angles.values.tolist() \
            and not forward in self.res2.angles.values.tolist() \
            and not reverse in self.res2.angles.values.tolist():
                zero_angles.append([self.q_id(at, q_atoms) for at in angle])

        return zero_angles

    def torsion_types(self):
        torsion_types = []
        return

    def zero_torsions(self, q_atoms):
        zero_torsions = []
        for torsion in self.hybrid.fep_torsions:
            forward = [x.upper() for x in torsion]
            reverse = forward[::-1]
            if not forward in self.res1.torsions.values.tolist() \
            and not reverse in self.res1.torsions.values.tolist() \
            and not forward in self.res2.torsions.values.tolist() \
            and not reverse in self.res2.torsions.values.tolist():
                zero_torsions.append([self.q_atnum(at, q_atoms) for at in torsion])
        return zero_torsions

    def q_atnum(self, atom, q_atoms):
        q_atnum = q_atoms[q_atoms.atom == atom].id.values[0]
        return q_atnum

    def q_id(self, atom, q_atoms):
        q_id = q_atoms[q_atoms.atom == atom].index.values[0]+1
        return q_id

    def write_fepfile(self, mutation, q_atoms, charges, atom_types, change_types, excluded_pairs, softcore, zero_angles):
        with open('tmp/fep1.fep', 'w') as fepfile:
            fepfile.write(f'! Hybrid topology FEP: {mutation}\n\n')
            fepfile.write(f'[FEP]\nstates 2\nsoftcore_use_max_potential on\n\n')
            fepfile.write(f'[atoms]\n')
            for i, atom in q_atoms.iterrows():
                fepfile.write(f'{i+1:>4} {atom.id:>6} ! {atom.atom}\n')
            fepfile.write(f'\n[change_charges]\n')
            for i, charge in charges.iterrows():
                fepfile.write(f'{i+1:>4} {charge.st1:>8} {charge.st2:>8}\n')
            fepfile.write(f'\n[atom_types]\n')
            for i, type in atom_types.iterrows():
                fepfile.write(f'{type.type:<4} {type.A1:>8} {type.B1:>6}  0.00  0.00 {type.A3:>8} {type.B2:>6} {type.mass:>6}\n')
            fepfile.write(f'\n[change_atoms]\n')
            for i, atom in change_types.iterrows():
                fepfile.write(f'{i+1:>4} {atom.st1:>4} {atom.st2:>4}\n')
            fepfile.write(f'\n[excluded_pairs]\n')
            for i, pair in excluded_pairs.iterrows():
                fepfile.write(f'{pair.at1:>6}{pair.at2:>6}  1  1\n')
            fepfile.write(f'\n[softcore]\n')
            for i, atom in softcore.iterrows():
                fepfile.write(f'{i+1:>4} {atom.st1:>4} {atom.st2:>4}\n')
            fepfile.write(f'\n[change_angles]\n')
            for i, angle in zero_angles.iterrows():
                fepfile.write(f'{angle.at1:>4} {angle.at2:>4} {angle.at3:>4}    0  0\n')
            fepfile.close()
        return

    def execute(self):
        # print(self.hybrid.hybrid_atoms)
        q_atoms = pd.DataFrame(self.get_top(), columns=['id', 'atom'])
        # print(q_atoms)
        charges = pd.DataFrame(self.charges(), columns=['st1', 'st2'])
        atom_types = pd.DataFrame(self.atom_types(), columns=['type', 'A1', 'A2', 'B1', 'A3', 'B2', 'mass'])
        change_types = pd.DataFrame(self.change_atoms(), columns=['st1', 'st2'])
        excluded_pairs = pd.DataFrame(self.excluded_pairs(q_atoms), columns=['at1', 'at2'])
        softcore = pd.DataFrame(self.softcore(), columns=['st1', 'st2'])
        zero_angles = pd.DataFrame(self.zero_angles(q_atoms), columns=['at1', 'at2', 'at3'])
        zero_torsions = pd.DataFrame(self.zero_torsions(q_atoms), columns=['at1', 'at2', 'at3', 'at4'])
        self.write_fepfile(self.mutation, q_atoms, charges, atom_types, change_types, excluded_pairs, softcore, zero_angles)
