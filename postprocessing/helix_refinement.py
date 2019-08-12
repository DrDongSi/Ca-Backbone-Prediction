"""Provides postprocessing step which refines the helices of the predicted
backbone

The idea behind the refinement is to adjust the nodes that are within a helix
in such way that they better represent the structure of an α-helix.

The helix refinement can be applied on a pdb file by calling the function
'apply_helix_refinement'.
"""

import numpy as np
from math import sin, cos
from scipy.optimize import minimize
from math import pi
from .pdb_reader_writer import PDB_Reader_Writer, Chain

__author__ = 'Jonas Pfab'


def update_paths(paths):
    paths['traces_refined'] = paths['output'] + 'traces_refined.pdb'


def execute(paths):
    """Coordinates the application of the helix refinement and writes new pdb
    file containing refined backbone structure"""
    reader_writer = PDB_Reader_Writer()
    chains = reader_writer.read_pdb(paths['traces_merged'])
    fit_helices(chains)
    reader_writer.write_pdb(chains, paths['traces_refined'])


def fit_helices(chains):
    """Fits nodes that belong to a helix to helix structure

    Parameters
    ----------
    chains: list
        List of 'Chain' objects containing information about nodes, sheets, and
        helices
    """
    for chain in chains:
        i = 0
        while i < len(chain.helices):
            helix_start, helix_end = chain.helices[i]
            helix_end = helix_end + 1
            original_len = len(chain.nodes)
            # Don't fit helices that consist of less than 10 atoms
            if helix_end - helix_start < 10:
                i += 1
                continue

            try:
                helix = Helix(interval_size=9, min_interval_size=3, r=2.11, c=1.149, gap=1.498, flatten=4)
                helix.fit(chain.nodes[helix_start:helix_end])

                chain.nodes[helix_start:helix_end] = helix.nodes
                diff = len(chain.nodes) - original_len

                update_sec_structure(chain.helices, helix_start, diff)
                update_sec_structure(chain.sheets, helix_start, diff)

                i += 1
            except CurvedScrewAxisError as error:
                split_helix_at_node(error.node, chain, i)


def update_sec_structure(sec_structure, node_start, diff):
    """Updates start and end indices of secondary structure information

    Used to update the indices when the length of a helix is changed due to the
    helix refinement.

    Parameters
    ----------
    sec_structure: list
        List of chains with secondary structure info for every chain

    node_start: int
        Index of node at which the indices need to be adjusted

    diff: int
        Difference by which the indices need to be updated
    """
    for node_range in sec_structure:
        if node_range[0] > node_start:
            node_range[0] += diff
        if node_range[1] > node_start:
            node_range[1] += diff


def split_helix_at_node(node, chain, i):
    """Splits helix from 'chain' at given index 'i' at closest node to 'node'

    Parameters
    ----------
    node: list
        Closest node in the helix to the node is where the helix is split.
        Since node is from the screw axis of the helix it cannot be the node
        itself

    chain: Chain
        Chain in which the helix is split

    i: int
        Index of helix is split
    """
    helix_start, helix_end = chain.helices[i]
    # Split helix into two helices at closest node
    closest_node, min_distance = -1, float('inf')
    for j in range(helix_start, helix_end):
        distance = get_distance(chain.nodes[j], node)
        if distance < min_distance:
            closest_node = j
            min_distance = distance

    if helix_start == closest_node:
        closest_node += 1

    del chain.helices[i]
    chain.helices.insert(i, [closest_node, helix_end])
    chain.helices.insert(i, [helix_start, closest_node])


class CurvedScrewAxisError(Exception):
    """Is raised when the screw axis of a helix exceeds a certain angle

    Parameters
    ----------
    node: list
        Node from screw axis at which the max angle is exceeded
    """

    def __init__(self, node):
        super(CurvedScrewAxisError, self).__init__('Screw axis is too curved at node: ' + str(node))
        self.node = node


class Helix:
    """Represents a helix along a screw axis

    Parameters
    ----------
    interval_size: int
        Number of consecutive nodes used to calculate a centroid

    min_interval_size: int
        Number of consecutive nodes required to calculate a centroid (The edges
        of the helix might be cleaner if this value is odd)

    r: float
        Defines the radius of the helix

    c: float
        Defines the vertical separation of the helix's loops

    gap: float
        Gap between atoms when traveling along the screw axis

    flatten: int
        Determines how much the screw axis of the helix is flattened

    Notes
    ----------
    The Helix class does not automatically represent an Alpha-Helix. Only with
    the appropriate parameters will it resemble one.
    """

    def __init__(self, interval_size, min_interval_size, r, c, gap, flatten):
        self.screw_axis = None
        self.nodes = []
        self.interval_size = interval_size
        self.min_interval_size = min_interval_size
        self.r = r
        self.c = c
        self.gap = gap
        self.flatten = flatten

    def fit(self, nodes):
        """Fits Helix onto given list of nodes

        Method automatically finds screw axis of the helix based on the nodes.
        Next, it draws the helix with different shift/rotation values and uses
        the one which is closest to the original list of nodes.

        Parameters
        ----------
        nodes: list
            List of coordinates of the predicted location of the CA atoms
        """
        self._set_screw_axis(nodes)

        node_max_angle = self.screw_axis.get_node_max_angle(80, interval_size=3)
        if node_max_angle is not None:
            raise CurvedScrewAxisError(node_max_angle)

        def evaluate_params(x):
            self._set_nodes(shift=x[0], rotation=x[1])

            return get_avg_offset(nodes, self.nodes)

        # Find best shift and rotation parameters
        # noinspection PyTypeChecker
        res = [minimize(evaluate_params, [-self.gap, -2 * pi]),
               minimize(evaluate_params, [0, 0]),
               minimize(evaluate_params, [self.gap, 2 * pi])]

        res.sort(key=lambda params: evaluate_params(params.x))

        self._set_nodes(shift=res[0].x[0], rotation=res[0].x[1])
        if 0.8 < get_avg_offset(nodes, self.nodes) <= 1.7:
            self._set_edges(nodes)
        else:
            self.nodes = nodes

    def _set_screw_axis(self, nodes):
        """Finds screw axis by calculating centroids for different intervals"""
        centroids = []
        interval = []
        for node in nodes:
            interval.append(node)

            # Make sure the centroid is always calculated by the last interval
            # size elements
            if len(interval) > self.interval_size:
                del interval[0]

            if len(interval) >= self.min_interval_size:
                centroids.append(get_centroid(interval))

        while len(interval) >= self.min_interval_size:
            del interval[0]
            centroids.append(get_centroid(interval))

        self.screw_axis = Curve(centroids, self.flatten)

    def _set_nodes(self, shift, rotation):
        """Set new nodes which conform to the helical structure

        The helix is build using the following formula:
        x(t) = r * sin((c * t) - shift + rotation)
        y(t) = r * cos((c * t) - shift + rotation)
        z(t) = t

        The vector going from the point at the screw axis at t to the point
        calculated by the formula is then mapped to the plane orthogonal to the
        direction of the screw axis at t.

        Parameters
        ----------
        shift: float
            Initial shift of the first CA atom when building the helix

        rotation: float
            Initial rotation of the helix
        """
        self.nodes = []
        t = shift
        while self.screw_axis(t) is not None:
            # The new axises to which the vector is mapped
            z_axis = -1 * normalize(self.screw_axis.get_vector(t))
            y_axis = -1 * normalize(np.cross(z_axis, np.array([1, 0, 0])))
            x_axis = -1 * normalize(np.cross(z_axis, y_axis))
            rotation_matrix = np.array([x_axis, y_axis, z_axis])

            original_vector = np.array([self.r * sin((self.c * t) - shift + rotation),
                                        self.r * cos((self.c * t) - shift + rotation),
                                        0])

            rotated_vector = np.dot(original_vector, rotation_matrix)

            self.nodes.append(np.add(self.screw_axis(t), rotated_vector))
            t += self.gap

    def _set_edges(self, original_nodes):
        """Sets the edge nodes which were ignored due to the min interval size"""
        edge_length = int(self.min_interval_size / 2 - 0.1)

        self.nodes = original_nodes[:edge_length] + self.nodes + original_nodes[-edge_length:]


class Curve:
    """A curve represented by interval lines

    Parameters
    ----------
    nodes: list
        List of nodes through which the curve should all go through

    flatten: int
        Determines how much the curve is flattened
    """

    def __init__(self, nodes, flatten):
        # Remove every flatten node but always keep first and last nodes
        self.nodes = nodes[:-1][::flatten] + [nodes[-1]]

    def __call__(self, t):
        """Returns the node from the curve at t"""
        for i in range(len(self.nodes) - 1):
            vector_to_next = get_vector(self.nodes[i], self.nodes[i + 1])
            distance_to_next = np.linalg.norm(vector_to_next)

            if distance_to_next > t:
                return self.nodes[i] + ((t / distance_to_next) * vector_to_next)
            else:
                t -= distance_to_next

        return None

    def get_vector(self, t):
        """Returns the direction vector of the curve at t"""
        for i in range(len(self.nodes) - 1):
            vector_to_next = get_vector(self.nodes[i], self.nodes[i + 1])
            distance_to_next = np.linalg.norm(vector_to_next)

            if distance_to_next > t:
                return vector_to_next
            else:
                t -= distance_to_next

        return None

    def get_node_max_angle(self, max_angle, interval_size):
        """Returns first node from where curve exceeds max angle

        Parameters
        ----------
        max_angle: float
            Maximum angle in degrees that curve can have within 'interval_size'
            nodes

        interval_size: int
            Number of nodes in interval for which the angle is compared to
            'max_angle'
        """
        current_interval = np.array([0] * interval_size)
        for i in range(1, len(self.nodes) - 1):
            vector_to_current = get_vector(self.nodes[i - 1], self.nodes[i])
            vector_to_next = get_vector(self.nodes[i], self.nodes[i + 1])

            current_interval[i % interval_size] = angle_between(vector_to_current, vector_to_next)

            if sum(current_interval) > max_angle:
                return self.nodes[max(i - int(interval_size / 2), 0)]

        return None


def get_avg_offset(nodes1, nodes2):
    """Calculates the average distance from nodes in nodes1 to the closest
     nodes in nodes2"""
    sum_offsets = 0
    for node1 in nodes1:
        closest_distance = -1
        for node2 in nodes2:
            distance = get_distance(node1, node2)
            if distance < closest_distance or closest_distance == -1:
                closest_distance = distance

        sum_offsets += closest_distance

    return sum_offsets / len(nodes1)


def get_distance(point1, point2):
    """Returns distance between point1 and point2"""
    return np.linalg.norm(point1 - point2)


def get_vector(node1, node2):
    """Creates vector which connects node1 with node2"""
    return np.array([node2[0] - node1[0], node2[1] - node1[1], node2[2] - node1[2]])


def normalize(v):
    """Normalizes vector to unit vector"""
    norm = np.linalg.norm(v)
    if norm == 0:
        return v

    return v / norm


def get_centroid(nodes):
    """Returns centroid node from given list of nodes"""
    n = len(nodes)
    sum_x, sum_y, sum_z = 0, 0, 0

    for node in nodes:
        sum_x += node[0]
        sum_y += node[1]
        sum_z += node[2]

    return np.array([sum_x / n, sum_y / n, sum_z / n])


def project_vector(vector1, vector2):
    """Projects vector1 onto vector2"""
    return (np.dot(vector1, vector2) / np.linalg.norm(vector2)) * (vector2 / np.linalg.norm(vector2))


def angle_between(v1, v2):
    """Returns the angle in radians between vectors 'v1' and 'v2'"""
    v1_u = normalize(v1)
    v2_u = normalize(v2)

    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) * (180 / pi)