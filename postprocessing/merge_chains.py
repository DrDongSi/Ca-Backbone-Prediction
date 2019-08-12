import numpy as np
from .pdb_reader_writer import PDB_Reader_Writer, Chain


def update_paths(paths):
    paths['traces_merged'] = paths['output'] + 'traces_merged.pdb'


def execute(paths):
    reader_writer = PDB_Reader_Writer()
    chains = [c for c in reader_writer.read_pdb(paths['traces']) if len(c.nodes) > 0]
    while merge_closest_chains(chains):
        pass

    chains = [c for c in chains if len(c.nodes) > 5]

    reader_writer.write_pdb(chains, paths['traces_merged'])


class PossibleConnection:

    def __init__(self, chain1, chain2, i1, i2):
        self.chain1 = chain1
        self.chain2 = chain2
        self.i1 = i1
        self.i2 = i2
        self.distance = get_distance(chain1.nodes[i1], chain2.nodes[i2])


def merge_closest_chains(chains):
    possible_connections = []
    for chain1 in chains:
        for chain2 in [c for c in chains if c != chain1]:
            last_i1, last_i2 = len(chain1.nodes) - 1, len(chain2.nodes) - 1
            possible_connections += [
                PossibleConnection(chain1, chain2, 0, 0),
                PossibleConnection(chain1, chain2, 0, last_i2),
                PossibleConnection(chain1, chain2, last_i1, 0),
                PossibleConnection(chain1, chain2, last_i1, last_i2)
            ]

    possible_connections.sort(key=lambda p: p.distance)
    possible_connections = list(filter(lambda p: p.distance > 3.2, possible_connections))
    if len(possible_connections) > 0 and possible_connections[0].distance < 10:
        merge_chains(chains,
                     possible_connections[0].chain1,
                     possible_connections[0].chain2,
                     possible_connections[0].i1,
                     possible_connections[0].i2)

        return True
    else:
        return False


def merge_chains(chains, chain1, chain2, at1, at2):
    if at1 == 0 and at2 == 0:
        chain1.helices = reverse_indices(chain2.helices, chain2.nodes) + add_offset(chain1.helices, len(chain2.nodes))
        chain1.sheets = reverse_indices(chain2.sheets, chain2.nodes) + add_offset(chain1.sheets, len(chain2.nodes))
        chain1.nodes = chain2.nodes[::-1] + chain1.nodes
        chains.remove(chain2)
    elif at1 == 0:
        chain2.helices += add_offset(chain1.helices, len(chain2.nodes))
        chain2.sheets += add_offset(chain1.sheets, len(chain2.nodes))
        chain2.nodes += chain1.nodes
        chains.remove(chain1)
    elif at2 == 0:
        chain1.helices += add_offset(chain2.helices, len(chain1.nodes))
        chain1.sheets += add_offset(chain2.sheets, len(chain1.nodes))
        chain1.nodes += chain2.nodes
        chains.remove(chain2)
    else:
        chain1.helices += add_offset(reverse_indices(chain2.helices, chain2.nodes), len(chain1.nodes))
        chain1.sheets += add_offset(reverse_indices(chain2.sheets, chain2.nodes), len(chain1.nodes))
        chains.remove(chain2)
        chain1.nodes += chain2.nodes[::-1]


def add_offset(sse, offset):
    return [[i + offset for i in i_s] for i_s in sse]


def reverse_indices(sse, nodes):
    return [[len(nodes) - i - 1 for i in i_s][::-1] for i_s in sse][::-1]


def get_distance(point1, point2):
    """Returns distance between point1 and point2"""
    return np.linalg.norm(point1 - point2)