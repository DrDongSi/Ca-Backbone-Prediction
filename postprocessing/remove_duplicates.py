"""Removes duplicate nodes (same coordinates)"""

from .pdb_reader_writer import PDB_Reader_Writer


def update_paths(paths):
    paths['duplicates_removed'] = paths['output'] + 'duplicates_removed.pdb'


def execute(paths):
    reader_writer = PDB_Reader_Writer()
    chains = reader_writer.read_pdb(paths['traces_refined'])
    remove_duplicates(chains)
    reader_writer.write_pdb(chains, paths['duplicates_removed'])


def remove_duplicates(chains):
    for i in range(len(chains)):
        chain = chains[i]
        for other_chain in (chains[0:i] + chains[i + 1:]):
            for j, k in [(0, 0), (0, -1), (-1, 0), (-1, -1)]:
                if not chain.nodes or not other_chain.nodes:
                    continue

                if is_equal(chain.nodes[j], other_chain.nodes[k]):
                    if len(chain.nodes) > len(other_chain.nodes):
                        del other_chain.nodes[k]
                    else:
                        del chain.nodes[j]


def is_equal(node1, node2):
    return node1[0] == node2[0] and node1[1] == node2[1] and node1[2] == node2[2]
