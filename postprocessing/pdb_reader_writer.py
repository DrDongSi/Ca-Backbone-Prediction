import numpy as np

########################################################################################################################
# Data IO methods
########################################################################################################################

class Chain:
    """Tracks nodes, sheets, and helices for a certain chain"""
    def __init__(self):
        self.nodes = []
        self.sheets = []
        self.helices = []

class PDB_Reader_Writer:
    def read_pdb(self, pdb_name):
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
                        data = self.__parse_node(line)
                        chains[-1].nodes.append(data)
                    elif line[:5] == 'HELIX':
                        data, i = self.__parse_helix(line, chains)
                        chains[i].helices.append(data)
                    elif line[:5] == 'SHEET':
                        data, i = self.__parse_sheet(line, chains)
                        chains[i].sheets.append(data)
                except ValueError as error:
                    print('Error parsing ' + line + str(error))

        return chains


    def write_pdb(self, chains, file_name):
        """Writes nodes and secondary structure info to pdb file with given name"""
        nodes_str = []
        helices_str = []
        sheets_str = []

        offset = 0
        for chain in chains:
            for i in range(len(chain.nodes)):
                nodes_str.append(self.__format_node(chain.nodes[i], 'A', offset + i + 1))
            nodes_str.append('TER\n')

            for helix_start, helix_end in chain.helices:
                helices_str.append(self.__format_helix_info('A', helix_start + offset, helix_end + offset))

            for sheet_start, sheet_end in chain.sheets:
                sheets_str.append(self.__format_sheet_info('A', sheet_start + offset, sheet_end + offset))

            offset += len(chain.nodes) + 1

        with open(file_name, 'w') as pdb_file:
            pdb_file.write(''.join(nodes_str))
            pdb_file.write(''.join(helices_str))
            pdb_file.write(''.join(sheets_str))


    def __parse_node(self, line):
        """Parses node data from given 'line'"""
        return np.array([float(line[30:38]),
                         float(line[38:46]),
                         float(line[46:54])])


    def __parse_helix(self, line, chains):
        """Parses helix data from given 'line' and calculates chain index 'i' from
        'chains'"""
        i = 0
        data = [int(line[21:25]) - 1, int(line[33:37]) - 1]
        for chain in chains:
            if data[1] <= len(chain.nodes):
                break

            data[0] -= len(chain.nodes) + 1
            data[1] -= len(chain.nodes) + 1
            i += 1

        return data, i


    def __parse_sheet(self, line, chains):
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


    def __format_node(self, node, chain, n):
        """Encodes node to str in PDB format"""
        return \
            'ATOM      1  CA  GLY ' + \
            chain + \
            str(n).rjust(4) + '    ' + \
            '{0:.3f}'.format(node[0]).rjust(8) + \
            '{0:.3f}'.format(node[1]).rjust(8) + \
            '{0:.3f}'.format(node[2]).rjust(8) + \
            '  1.00  0.00           C  \n'


    def __format_helix_info(self, chain, node_from, node_to):
        """Encodes helix info to str in PDB format"""
        return \
            'HELIX    1   1 GLY ' + \
            chain + ' ' + \
            str(node_from + 1).rjust(4) + '  GLY ' + \
            chain + ' ' + \
            str(node_to + 1).rjust(4) + '  1\n'


    def __format_sheet_info(self, chain, node_from, node_to):
        """Encodes sheet info to str in PDB format"""
        return \
            'SHEET    1   A 6 GLY ' + \
            chain + \
            str(node_from + 1).rjust(4) + '  GLY ' + \
            chain + \
            str(node_to + 1).rjust(4) + '  0\n'



