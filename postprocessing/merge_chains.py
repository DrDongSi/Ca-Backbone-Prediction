import numpy as np
from collections import defaultdict


def update_paths(paths):
    paths['fragments_merged'] = paths['output'] + 'fragments_merged.pdb'


def execute(paths):
    chains = [c for c in read_pdb(paths['traces_refined']) if len(c.nodes) > 0]
    #while merge_closest_chains(chains):
        #pass
    merge_closest_chains(chains)

    write_pdb(chains, paths['fragments_merged'])


class EndPoint:
    def __init__(self, chains, chain1, chain2, closest_nodes1, closest_nodes2, distance):
        self.chains = chains
        self.chain1 = chain1
        self.chain2 = chain2
        self.closest_nodes0 = closest_nodes1
        self.closest_nodes1 = closest_nodes2
        self.distance = distance
    def __str__(self):
        return str(self.distance)

def merge_closest_chains(chains):


    matches = defaultdict(list)
    for chain1 in chains:
        for chain2 in [c for c in chains if c != chain1]:
            last_i1, last_i2 = len(chain1.nodes) - 1, len(chain2.nodes) - 1
            closest_nodes = sorted([
                (0, 0, get_distance(chain1.nodes[0], chain2.nodes[0])),
                (0, last_i2, get_distance(chain1.nodes[0], chain2.nodes[last_i2])),
                (last_i1, 0, get_distance(chain1.nodes[last_i1], chain2.nodes[0])),
                (last_i1, last_i2, get_distance(chain1.nodes[last_i1], chain2.nodes[last_i2]))
            ], key=lambda d: d[2])[0]

            if closest_nodes[2] < 10:
                print(closest_nodes)
                distance = closest_nodes[2]
                matches[chain1].append(EndPoint(chains, chain1, chain2, closest_nodes[0], closest_nodes[1], distance))
                #merge_chains(chains, chain1, chain2, closest_nodes[0], closest_nodes[1])
                #return True

    for key in matches.keys():
        matches[key] = sorted(matches[key], key=lambda d: d.distance)
        matches[key] = list(filter(lambda d: d.distance > 0, matches[key]))
        print('--->: ')
        for item in matches[key]:
            print(item)

    for key, item in matches.items():
        if len(item) > 0:
            merge_chains(item[0].chains, item[0].chain1, item[0].chain2, item[0].closest_nodes0, item[0].closest_nodes1)

    #return False


def merge_chains(chains, chain1, chain2, at1, at2):
    if at1 == 0 and at2 == 0:
        chain1.nodes = chain2.nodes[::-1] + chain1.nodes
        chain1.helices = reverse_indices(chain2.helices, chain2.nodes) + add_offset(chain1.helices, len(chain1.nodes))
        chain1.sheets = reverse_indices(chain2.sheets, chain2.nodes) + add_offset(chain1.sheets, len(chain1.nodes))
        #chains.remove(chain2)
    elif at1 == 0:
        chain2.nodes += chain1.nodes
        chain2.helices += add_offset(chain1.helices, len(chain2.nodes))
        chain2.sheets += add_offset(chain1.sheets, len(chain2.nodes))
        #chains.remove(chain1)
    elif at2 == 0:
        chain1.nodes += chain2.nodes
        chain1.helices += add_offset(chain2.helices, len(chain2.nodes))
        chain1.sheets += add_offset(chain2.sheets, len(chain2.nodes))
        #chains.remove(chain2)
    else:
        chain1.nodes += chain2.nodes[::-1]
        chain1.helices += add_offset(reverse_indices(chain2.helices, chain2.nodes), len(chain1.nodes))
        chain1.sheets += add_offset(reverse_indices(chain2.sheets, chain2.nodes), len(chain1.nodes))
        #chains.remove(chain2)


def add_offset(sse, offset):
    return [[i + offset for i in i_s] for i_s in sse]


def reverse_indices(sse, nodes):
    return [[len(nodes) - i for i in i_s] for i_s in sse]


def get_distance(point1, point2):
    """Returns distance between point1 and point2"""
    return np.linalg.norm(point1 - point2)


########################################################################################################################
# Data IO methods
########################################################################################################################


class Chain:
    """Tracks nodes, sheets, and helices for a certain chain"""
    def __init__(self):
        self.nodes = []
        self.sheets = []
        self.helices = []


def read_pdb(pdb_name):
    """Reads pdb file with given name

    Parameters
    ----------
    pdb_name: str
        Name of pdb file

    Returns
    ----------
    chains: list
        List of 'Chain' objects which contains nodes and sheets information

    Notes
    ----------
    Information about helices is not parsed from the pbd file but rather has to
    be parsed from the helix mrc file
    """
    chains = [Chain()]
    with open(pdb_name) as pdb_file:
        for line in pdb_file:
            if 'TER' in line:
                chains.append(Chain())
            try:
                if line[:4] == 'ATOM':
                    data = parse_node(line)
                    chains[-1].nodes.append(data)
                elif line[:5] == 'HELIX':
                    data, i = parse_helix(line, chains)
                    chains[i].helices.append(data)
                elif line[:5] == 'SHEET':
                    data, i = parse_sheet(line, chains)
                    chains[i].sheets.append(data)
            except ValueError as error:
                print('Error parsing ' + line + str(error))

    return chains


def write_pdb(chains, file_name):
    """Writes nodes and secondary structure info to pdb file with given name"""
    nodes_str = []
    helices_str = []
    sheets_str = []

    offset = 0
    for chain in chains:
        for i in range(len(chain.nodes)):
            nodes_str.append(format_node(chain.nodes[i], 'A', offset + i + 1))
        nodes_str.append('TER\n')

        for helix_start, helix_end in chain.helices:
            helices_str.append(format_helix_info('A', helix_start + offset, helix_end + offset))

        for sheet_start, sheet_end in chain.sheets:
            sheets_str.append(format_sheet_info('A', sheet_start + offset, sheet_end + offset))

        offset += len(chain.nodes) + 1

    with open(file_name, 'w') as pdb_file:
        pdb_file.write(''.join(nodes_str))
        pdb_file.write(''.join(helices_str))
        pdb_file.write(''.join(sheets_str))


def parse_node(line):
    """Parses node data from given 'line'"""
    return np.array([float(line[31:38]),
                     float(line[39:46]),
                     float(line[47:54])])


def parse_helix(line, chains):
    """Parses helix data from given 'line' and calculates chain index 'i' from
    'chains'"""
    i = 0
    data = [int(line[21:25]) - 1, int(line[33:37])]
    for chain in chains:
        if data[1] <= len(chain.nodes):
            break

        data[0] -= len(chain.nodes) + 1
        data[1] -= len(chain.nodes) + 1
        i += 1

    return data, i


def parse_sheet(line, chains):
    """Parses sheet data from given 'line' and calculates chain index 'i' from
    'chains'"""
    i = 0
    data = [int(line[22:26]) - 1, int(line[33:37]) - 1]
    for chain in chains:
        if data[1] < len(chain.nodes):
            break

        data[0] -= len(chain.nodes) + 1
        data[1] -= len(chain.nodes) + 1
        i += 1

    return data, i


def format_node(node, chain, n):
    """Encodes node to str in PDB format"""
    return \
        'ATOM      1  CA  GLY ' + \
        chain + \
        str(n).rjust(4) + '    ' + \
        '{0:.3f}'.format(node[0]).rjust(8) + \
        '{0:.3f}'.format(node[1]).rjust(8) + \
        '{0:.3f}'.format(node[2]).rjust(8) + \
        '  1.00  0.00           C  \n'


def format_helix_info(chain, node_from, node_to):
    """Encodes helix info to str in PDB format"""
    return \
        'HELIX    1   1 GLY ' + \
        chain + ' ' + \
        str(node_from + 1).rjust(4) + '  GLY ' + \
        chain + ' ' + \
        str(node_to).rjust(4) + '  1\n'


def format_sheet_info(chain, node_from, node_to):
    """Encodes sheet info to str in PDB format"""
    return \
        'SHEET    1   A 6 GLY ' + \
        chain + \
        str(node_from + 1).rjust(4) + '  GLY ' + \
        chain + \
        str(node_to + 1).rjust(4) + '  0\n'
