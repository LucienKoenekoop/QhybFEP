import networkx as nx

class Graph(object):

    def __init__(self, hybrid):
        graph = nx.Graph(nx.from_pandas_edgelist(hybrid.hybrid_bonds, 'atom1', 'atom2'))

        angles = []
        fep_angles = []
        torsions = []
        fep_torsions = []
        for i, atom_i in enumerate(list(graph.nodes)):
            for atom_j in list(graph.nodes)[i+1:]:
                path = nx.shortest_path(graph, atom_i, atom_j)
                if len(path) == 3:
                    if not all([i.islower() for i in path]) and not all([i.isupper() for i in path]):
                        fep_angles.append(path.copy())
                if len(path) == 4:
                    if not all([i.islower() for i in path]) and not all([i.isupper() for i in path]):
                        fep_torsions.append(path.copy())
                for x, atom in enumerate(path):
                    path[x] = hybrid.atom2type(atom)
                reverse = path[::-1]
                if len(path) == 3:
                    if not path in angles and not reverse in angles:
                        angles.append(path)
                elif len(path) == 4:
                    if not path in torsions and not reverse in torsions:
                        torsions.append(path)
                else:
                    continue

        hybrid.hybrid_angles = angles
        hybrid.hybrid_torsions = torsions
        hybrid.fep_angles = fep_angles
        hybrid.fep_torsions = fep_torsions
        return