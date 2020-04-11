"""This post-processing step will create a pdb file with the backbone trace from
the predictions of the CNN

The two main parts of this step is the 'confidence_walk' method, which takes the
output from the CNN and walks the confidence regions of Ca and Backbone
prediction to produce a series of traces that have a high-confidence of being
the actual backbone in the protein, and the 'Graph' class through which the
traces are refined.
"""

import mrcfile
import numpy as np
from copy import deepcopy
import math
from collections import deque
from .pdb_reader_writer import PDB_Reader_Writer

__author__ = 'Spencer Moritz'


def update_paths(paths):
    paths['first_confidence_walk'] = paths['output'] + 'first_confidence_walk.pdb'
    paths['second_confidence_walk'] = paths['output'] + 'second_confidence_walk.pdb'
    paths['refined_backbone'] = paths['output'] + 'refined_backbone.mrc'
    paths['final_ca_prediction'] = paths['output'] + 'final_ca_prediction.pdb'
    paths['traces'] = paths['output'] + 'traces.pdb'


def execute(paths):
    # Open up the required files.
    normalized_map = mrcfile.open(paths['normalized_map'], mode='r')
    ca_image = mrcfile.open(paths['ca_confidence'], mode='r').data
    helix_image = mrcfile.open(paths['helix_confidence'], mode='r').data
    sheet_image = mrcfile.open(paths['sheet_confidence'], mode='r').data
    backbone_image = mrcfile.open(paths['backbone_confidence'], mode='r').data
    origin = normalized_map.header.origin.item(0)

    # This is the major Post-Processing step where the full backbone trace is
    # built. This may take a few minutes to execute.
    confidence_walk(deepcopy(ca_image), origin, backbone_image, paths['first_confidence_walk'])

    # Build graph and then clean it to improve final output.
    graph = make_graph(paths['first_confidence_walk'])
    graph.edge_check()
    graph.remove_side_chains()
    graph.remove_loops(ca_image, origin)
    graph.remove_single_ends()
    graph.remove_empty_nodes()
    graph.remove_side_chains()
    graph.remove_loops(ca_image, origin)
    graph.remove_single_ends()
    graph.remove_empty_nodes()

    new_backbone = graph.refine_backbone(backbone_image, origin)
    with mrcfile.new(paths['refined_backbone'], overwrite=True) as mrc:
        mrc.set_data(new_backbone)
        mrc.header.origin = normalized_map.header.origin.item(0)
        mrc.update_header_stats()
        mrc.close()

    confidence_walk(deepcopy(ca_image), origin, new_backbone, paths['second_confidence_walk'])

    # Build graph and then clean it to improve final output.
    graph = make_graph(paths['second_confidence_walk'])
    graph.edge_check()
    graph.remove_side_chains()
    graph.remove_loops(ca_image, origin)
    graph.remove_single_ends()
    graph.remove_empty_nodes()
    graph.remove_side_chains()
    graph.remove_loops(ca_image, origin)
    graph.remove_single_ends()
    graph.remove_tail_loops()
    graph.remove_empty_nodes()

    # Now print the final graphs to output
    graph.print_graph(paths['final_ca_prediction'])
    graph.print_traces(sheet_image, helix_image, origin, paths['traces'])


def confidence_walk(prediction_image, offset, backbone_image, output_file):
    """ It produces a series of traces representative of the protein's backbone
    structure

    This is the main function used for the path-walking.

    prediction_image: the Ca-confidence 3D image.
    offset: A 3x1 vectors containing the x, y, z offset of the original .PDB image.
    pdbId: ID of the given protein, used to prefix the file name when printing.
    ss_image: A 3D image containing the a-helix confidence image.
    backbone_image: A 3D image containing the backbone confidence prediction of the protein.
    The primary loop of this function path-walks the prediction image along the backbone
    prediction. It finds areas of high confidence for Ca atoms and places an atom at each
    location. Once each trace has been founds, this function connects each trace and
    then prints the final graph to a file for further processing.
    """
    num_ca_edges_hash = np.zeros((np.shape(prediction_image)))
    untouched_prediction = deepcopy(prediction_image)
    set_of_ca_sets = list()
    for index in range(2436111 + 1):
        # Find and update for the high-confident location
        location = find_highest_confidence_ca(prediction_image, set_of_ca_sets)
        if location is None:
            break
        update_confidence_image(prediction_image, num_ca_edges_hash, location)

        # Find and update for the neighbor
        neighbor = find_nearest_ca(location, backbone_image, set_of_ca_sets, untouched_prediction)
        if neighbor is not None:
            update_confidence_image(prediction_image, num_ca_edges_hash, neighbor)
            update_ca_sets(set_of_ca_sets, location, neighbor)
            print('Placed Edge: ' + str(index + 1))

        if index % 10 == 0:
            print_ca_sets(set_of_ca_sets, offset, output_file)

    massage_ends(set_of_ca_sets)
    overlay_cas(set_of_ca_sets)

    print_ca_sets(set_of_ca_sets, offset, output_file)  # Final print


def build_graph(pdb_file, ca_image, origin, remove_tail_loops):
    """Build graph and then clean it to improve final output"""
    graph = make_graph(pdb_file)
    graph.edge_check()
    graph.remove_side_chains()
    graph.remove_loops(ca_image, origin)
    graph.remove_single_ends()
    graph.remove_empty_nodes()
    graph.remove_side_chains()
    graph.remove_loops(ca_image, origin)
    graph.remove_single_ends()
    if remove_tail_loops:
        graph.remove_tail_loops()
    graph.remove_empty_nodes()


def find_highest_confidence_ca(remaining_image, set_of_ca_sets):
    """This function finds the next location in the entire protein image to
    keep path-walking from

    If there is an active trace in the current path, then this function will return that point.
    If all current traces have been terminated, then this function will return the point
    in the prediction with the highest intensity to start a new path-walk.
    """
    max_location = None
    max_value = 0
    for ca_set in set_of_ca_sets:
        for index in range(len(ca_set)):
            if index == 0 or index == len(ca_set) - 1:
                value = remaining_image[ca_set[index][0], ca_set[index][1], ca_set[index][2]]
                if value > max_value:
                    max_value = value
                    max_location = ca_set[index]
    if max_value <= 0:
        # Only do this if we cannot find a good fit on the current back-chain
        box_size = np.shape(remaining_image)
        argmax = np.argmax(remaining_image)
        location = [0] * 3
        location[0] = int(argmax / (box_size[1] * box_size[2]))
        location[1] = int((argmax - location[0] * box_size[1] * box_size[2]) / box_size[2])
        location[2] = int(argmax - location[0] * box_size[1] * box_size[2] - location[1] * box_size[2])
        max_location = location
        # Quitting Condition
        if remaining_image[max_location[0], max_location[1], max_location[2]] <= 8:
            return None

    return max_location


def update_confidence_image(prediction_image, num_ca_edges_hash, location):
    """Updates the prediction image by zeroing out space that already has a
    placed Ca-atom"""
    num_ca_edges_hash[location[0]][location[1]][location[2]] += 1
    level = num_ca_edges_hash[location[0]][location[1]][location[2]]
    zero_out_sphere(prediction_image, location, level)


def zero_out_sphere(remaining_image, location, level):
    """This function zeros out a sphere of 3A in the remaining image

    This is necessary to prevent the path-walking algorithm from using repeated
    locations in the image. Essentially, this is a way of making a given space
    in the image as 'visited'.
    """
    box_size = np.shape(remaining_image)
    sphere_radius = 3
    for x in range(-sphere_radius + location[0], sphere_radius + location[0]):
        for y in range(-sphere_radius + location[1], sphere_radius + location[1]):
            for z in range(-sphere_radius + location[2], sphere_radius + location[2]):
                if (0 <= x < box_size[0] and 0 <= y < box_size[1] and 0 <= z < box_size[2] and
                        distance(location, (x, y,z)) < sphere_radius ** 2):
                    if level == 2 or not (x == location[0] and y == location[1] and z == location[2]):
                        remaining_image[x][y][z] = 0


def distance(x, y):
    """Calculates Euclidean distance between two points"""
    x = np.array(x)
    y = np.array(y)
    delta = x - y
    d2 = (delta * delta).sum()
    return d2 ** 0.5


def find_nearest_ca(previous_location, backbone_image, set_of_ca_sets, untouched_prediction):
    """This method finds the best possible neighboring Ca-atom

    The first part of the function generates all locations within the
    appropriate distance that have non-zero confidence values. It then
    calculates if these values are connected to the current Ca-location using a
    BFS. Each potential Ca-location is then given a score which is mostly a
    function of Ca-prediction cylindrical density. The voxel with the highest
    score is assigned as the next neighbor. If there is no suitable neighbor,
    None is returned and the trace is terminated.
    """
    invalid_ca_spots = find_my_neighbors(set_of_ca_sets, previous_location)
    invalid_ca_spots.append(previous_location)
    box_size = np.shape(untouched_prediction)
    sphere_radius = 5
    max_ca_distance = 4.5
    min_ca_distance = 3
    possible_set = list()
    # Loop through neighboring voxels looking for locations that might be
    # suitable.
    for x in range(-sphere_radius + previous_location[0], sphere_radius + previous_location[0]):
        for y in range(-sphere_radius + previous_location[1], sphere_radius + previous_location[1]):
            for z in range(-sphere_radius + previous_location[2], sphere_radius + previous_location[2]):
                if (0 <= x < box_size[0] and 0 <= y < box_size[1] and 0 <= z < box_size[2] and
                        max_ca_distance >= distance(previous_location, np.array((x, y, z))) >= min_ca_distance and
                        not already_placed([x, y, z], invalid_ca_spots) and untouched_prediction[x][y][z] > 0):
                    possible_set.append([x, y, z])
    if len(possible_set) == 0:
        return None
    bfs_distance_image = distance_between_bfs(previous_location, backbone_image)

    dictionary = {}
    for coordinate in possible_set:
        voids, density = cylindrical_density_fast(coordinate, previous_location, backbone_image, untouched_prediction)
        angle = find_angle2(set_of_ca_sets, previous_location, coordinate)
        score = density
        if (bfs_distance_image[coordinate[0], coordinate[1], coordinate[2]] < 100) and angle > 70:
            dictionary[score] = [coordinate[0], coordinate[1], coordinate[2]]
    key_list = dictionary.keys()
    sorted_key_list = list(sorted(key_list, reverse=True))
    if len(sorted_key_list) == 0:
        return None
    if sorted_key_list[0] <= 0:
        return None

    return dictionary.get(sorted_key_list[0])


def find_my_neighbors(set_of_ca_sets, location):
    """Finds the nearest 5 neighbors from a trace in each direction

    This is useful for preventing the path-walking method from looping on itself
    """
    neighbor_list = list()
    for ca_set in set_of_ca_sets:
        for index in range(len(ca_set)):
            if np.array_equal(ca_set[index], location):
                if index + 1 < len(ca_set):
                    neighbor_list.append(ca_set[index + 1])
                if index + 2 < len(ca_set):
                    neighbor_list.append(ca_set[index + 2])
                if index + 3 < len(ca_set):
                    neighbor_list.append(ca_set[index + 3])
                if index + 4 < len(ca_set):
                    neighbor_list.append(ca_set[index + 4])
                if index + 5 < len(ca_set):
                    neighbor_list.append(ca_set[index + 5])
                if index > 0:
                    neighbor_list.append(ca_set[index - 1])
                if index > 1:
                    neighbor_list.append(ca_set[index - 2])
                if index > 2:
                    neighbor_list.append(ca_set[index - 3])
                if index > 3:
                    neighbor_list.append(ca_set[index - 4])
                if index > 4:
                    neighbor_list.append(ca_set[index - 5])

    return neighbor_list


def already_placed(location, invalid_ca_spots):
    """Returns true if the given location is within 3A of invalid spots in the
    prediction image"""
    sphere_radius = 3
    for ca in invalid_ca_spots:
        if distance(location, ca) < sphere_radius:
            return True
    return False


def distance_between_bfs(position, input_image):
    """Third function uses a breadth first search to find all voxels that are
    within 8 connected spaces of a given coordinate

    This is used to determine if a neighboring Ca is connected through the the
    backbone structure.
    """
    visited = np.zeros((np.shape(input_image)))
    distance_image = np.full((np.shape(input_image)), 100)
    position_queue = deque()
    distance_queue = deque()
    position_queue.append(position)
    distance_queue.append(0)
    visited[position[0]][position[1]][position[2]] = 1
    while len(distance_queue) > 0 and distance_queue[0] < 8:
        cur_position = position_queue.popleft()
        cur_distance = distance_queue.popleft()
        offsets = [[0, 0, -1], [0, 0, 1], [0, -1, 0], [0, 1, 0], [-1, 0, 0], [1, 0, 0]]
        for index in range(len(offsets)):
            x = cur_position[0] + offsets[index][0]
            y = cur_position[1] + offsets[index][1]
            z = cur_position[2] + offsets[index][2]
            if x < len(input_image) and y < len(input_image[0]) and z < len(input_image[0][0]):
                if input_image[x][y][z] > 0 and visited[x][y][z] == 0:
                    position_queue.append([x, y, z])
                    distance_queue.append(cur_distance + 1)
                    distance_image[x][y][z] = cur_distance + 1
                    visited[x][y][z] = 1
    return distance_image


def cylindrical_density(ca_1, ca_2, input_image, untouched_prediction):
    """Calculates the 1A-radius cylindrical density between two Ca atoms

    This function returns two values. The first is the average number of voids
    (zero valued voxels) within the cylinders. The second is the density of the
    entire cylinder.
    """
    box_size = np.shape(input_image)
    steps = 10
    total_density = 0
    number_of_voids = 0
    number_of_points = 0
    midpoints = list()
    x_step = (ca_2[0] - ca_1[0]) / steps
    y_step = (ca_2[1] - ca_1[1]) / steps
    z_step = (ca_2[2] - ca_1[2]) / steps
    for index in range(steps + 1):
        midpoints.append([ca_1[0] + x_step * index, ca_1[1] + y_step * index, ca_1[2] + z_step * index])
    for z in range(int(ca_1[2]) - 4, int(ca_1[2]) + 5):  # 4 is kinda arbitrary here
        for y in range(int(ca_1[1]) - 4, int(ca_1[1]) + 5):
            for x in range(int(ca_1[0]) - 4, int(ca_1[0]) + 5):
                placed = False
                for index in range(len(midpoints)):
                    if (0 <= z < box_size[2] and box_size[1] > y >= 0 <= x < box_size[0] and
                            distance(midpoints[index], (x, y, z)) < 1
                            and not placed):
                        placed = True
                        total_density += untouched_prediction[x][y][z]
                        number_of_points += 1
                        if input_image[x][y][z] <= 0:
                            number_of_voids += 1
    return (number_of_voids * 15) / number_of_points, total_density / number_of_points


def cylindrical_density_fast(ca_1, ca_2, input_image, untouched_prediction):
    """Calculates the 1A-radius cylindrical density between two Ca atoms

    This function returns two values. The first is the average number of voids
    (zero valued voxels) within the cylinders. The second is the density of the
    entire cylinder.
    """
    box_size = np.shape(input_image)
    steps = 10
    frac = np.arange(steps + 1).reshape(steps + 1, 1) / 10
    midpoints = ca_1 * frac + ca_2 * (1 - frac)
    midf = np.floor(midpoints).astype(np.int)
    total_density = 0
    number_of_voids = 0
    number_of_points = 0
    voxels = set()
    for x in range(8):
        offset = np.array([[x & 1, (x & 2) // 2, (x & 4) // 4]])
        points = midf + offset
        # print('offset', offset)
        # print('points', points)
        assert points.min() >= 0
        assert (points.max() < np.array([box_size])).all()
        deltas = midpoints - points
        for point, delta in zip(points, deltas):
            # print('XXX', point, delta, (delta ** 2).sum() < 1)
            if (delta ** 2).sum() < 1:
                voxels.add(tuple(point.tolist()))

    voxel_list = np.array(list(voxels))

    for x, y, z in voxel_list:
        total_density += untouched_prediction[x, y, z]
        number_of_points += 1
        if input_image[x, y, z] <= 0:
            number_of_voids += 1

    return (number_of_voids * 15) / number_of_points, total_density / number_of_points


def find_angle2(set_of_ca_sets, old_location, new_location):
    """Another function to calculate the angle between three Ca-atoms"""
    for ca_set in set_of_ca_sets:
        for index in range(len(ca_set)):
            if np.array_equal(ca_set[index], old_location):
                if index + 1 < len(ca_set):
                    a = np.array([ca_set[index + 1][0], ca_set[index + 1][1], ca_set[index + 1][2]])
                elif index > 0:
                    a = np.array([ca_set[index - 1][0], ca_set[index - 1][1], ca_set[index - 1][2]])
                else:
                    return 180
                b = np.array([old_location[0], old_location[1], old_location[2]])
                c = np.array([new_location[0], new_location[1], new_location[2]])
                ba = a - b
                bc = c - b

                cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
                angle = np.arccos(cosine_angle)

                return np.degrees(angle)
    return 180


def update_ca_sets(set_of_ca_sets, location, neighbor):
    location_ca_set = find_ca_in_set_of_traces(location, set_of_ca_sets)
    if location_ca_set is None:
        location_ca_set = list()
        location_ca_set.append(location)
        location_ca_set.append(neighbor)
        set_of_ca_sets.append(location_ca_set)
    else:
        if location_ca_set[0] == location:
            location_ca_set.insert(0, neighbor)
        else:
            location_ca_set.append(neighbor)


def find_ca_in_set_of_traces(cur_ca, set_of_traces):
    """Finds the trace that contains a given Ca coordinate. If no trace
    contains this coordinate, returns None."""
    for ca_set in set_of_traces:
        for ca in ca_set:
            if np.array_equal(ca, cur_ca):
                return ca_set
    return None


def print_ca_sets(set_of_ca_sets, offset, file_name):
    """Prints the final set of Ca-traces to a new file that will be later read
    in the graph post processing step. Each trace is assigned its own chain
    letter."""
    confidence_walk_pdb = open(file_name, "w")
    counter = 0
    for index in range(len(set_of_ca_sets)):
        next_set = set_of_ca_sets[index]
        counter += 1
        for ca in range(len(next_set)):
            PDB_Reader_Writer.write_single_pdb(file=confidence_walk_pdb, type='ATOM', chain='A', node=np.array([(next_set[ca][2] + offset[0]),(next_set[ca][1] + offset[1]),(next_set[ca][0] + offset[2])]), seqnum=counter)
            counter += 1
    confidence_walk_pdb.close()


def massage_ends(set_of_ca_sets):
    """This method is used to connect disjoint traces at specific nodes

    Essentially, this method created a graph of traces from a set of disjoint
    graphs. This is accomplished by reassigning the ends og traces to lie on
    top of nodes within a neighboring trace."""
    for ca_set in set_of_ca_sets:
        switch_start_ca = None
        switch_start_distance = 1000
        switch_end_ca = None
        switch_end_distance = 1000
        ca_start = ca_set[0]
        ca_end = ca_set[len(ca_set) - 1]
        for next_ca_set in set_of_ca_sets:
            if next_ca_set != ca_set:
                for ca in next_ca_set:
                    start_dist = distance(ca, ca_start)
                    end_dist = distance(ca, ca_end)
                    if start_dist < 3 and start_dist < switch_start_distance:
                        switch_start_ca = deepcopy(ca)
                        switch_start_distance = start_dist
                    if end_dist < 3 and end_dist < switch_end_distance:
                        switch_end_ca = deepcopy(ca)
                        switch_end_distance = end_dist
        if switch_start_ca is not None:
            ca_set[0] = switch_start_ca
        if switch_end_ca is not None:
            ca_set[len(ca_set) - 1] = switch_end_ca


def overlay_cas(set_of_ca_sets):
    for ca_set in set_of_ca_sets:
        for index, ca in enumerate(ca_set):
            for next_ca_set in set_of_ca_sets:
                for next_index, next_ca in enumerate(next_ca_set):
                    if ca != next_ca:
                        if distance(ca, next_ca) < 3:
                            print('Before: ' + str(next_ca_set[next_index]))
                            next_ca_set[next_index] = ca_set[index]
                            print('After: ' + str(next_ca_set[next_index]))


class Node:
    """Simple Node class defining the location of a Ca atom along with which Ca
    atoms may be connected to it. There may be any number of connected Ca atoms."""

    def __init__(self, location):
        self.location = location
        self.edges = list()

    def get_location(self):
        return self.location

    def add_edge(self, edge):
        self.edges.append(edge)

    def get_num_edges(self):
        return len(self.edges)

    def get_edges(self):
        return self.edges


class Graph:
    """Graph class that defines the full backbone-Ca structure of the protein

    This class contains functions used to clean and improve the predicted
    backbone structure as a final step of post-processing"""

    def __init__(self):
        self.nodes = list()

    def add_node(self, node):
        self.nodes.append(node)

    def contains_location(self, location):
        for node in self.nodes:
            if node.get_location() == location:
                return True
        return False

    def get_node(self, location):
        for node in self.nodes:
            if node.get_location() == location:
                return node

        raise ValueError('No node found with given location')

    def edge_check(self):
        for node in self.nodes:
            index = 0
            while index < node.get_num_edges():
                if node.get_edges()[index] == node.get_location():
                    node.get_edges().remove(node.get_edges()[index])
                    index -= 1
                index += 1
            for edge in node.get_edges():
                if edge == node.get_location():
                    print("Edge Failure!")

    def get_end_nodes(self):
        end_nodes = list()
        for node in self.nodes:
            if node.get_num_edges() <= 1:
                end_nodes.append(node)
        return end_nodes

    def remove_single_ends(self):
        """This method removes all nodes/edges that are connected to a node
        that have 3 or more edges. It also removes Nodes that are connected
        to other end Nodes (a pair of disconnected Ca atoms)."""
        for node in self.nodes:
            if node.get_num_edges() == 1: # This is a single connection
                neighbor_node = self.get_node(node.get_edges()[0])
                if neighbor_node.get_num_edges() >= 3 or neighbor_node.get_num_edges() == 1:
                    neighbor_node.edges.remove(node.get_location())
                    node.edges.remove(neighbor_node.get_location())
                else:
                    for node_2 in neighbor_node.get_edges():
                        next_node = self.get_node(node_2)
                        if next_node.get_num_edges() >= 3:
                            next_node.edges.remove(neighbor_node.get_location())
                            neighbor_node.edges.remove(next_node.get_location())
                            neighbor_node.edges.remove(node.get_location())
                            node.edges.remove(neighbor_node.get_location())

    def walk_until_trinary(self, location, visited):
        """A recursive method used to find the end of a trace

        An end of a trace is defined as the point where the trace ends or when
        the trace hits a Ca atom with 3 or more edges"""
        node = self.get_node(location)
        visited.append(location)
        if node.get_num_edges() == 2:
            for edge in node.get_edges():
                if edge not in visited:
                    visited = self.walk_until_trinary(edge, deepcopy(visited))
        return visited

    def walk_graph(self, location, visited, depth):
        """A recursive method similar to walk_until_trinary that finds the
        depth of trace. Also stops when it reaches a node that has 3 edges
        or a terminating node."""
        node = self.get_node(location)
        if node.get_num_edges() <= 2:
            visited.append(location)
            for edge in node.get_edges():
                if edge not in visited:
                    depth = self.walk_graph(edge, deepcopy(visited), depth + 1)
        return depth

    def remove_single_loop(self, ca_list):
        trinary_node = self.get_node(ca_list[0])
        trinary_node.edges.remove(ca_list[1])
        trinary_node.edges.remove(ca_list[len(ca_list) - 1])
        for index, ca in enumerate(ca_list):
            if index > 0:
                to_be_removed = self.get_node(ca_list[index])
                self.nodes.remove(to_be_removed)

    def remove_tail_loops(self):
        for node in self.nodes:
            if node.get_num_edges() >= 3: # This is a trinary connection
                walk_lists = list()
                visited = list()
                visited.append(node.get_location())
                for edge in node.get_edges():
                    walk_lists.append(self.walk_until_trinary(edge, deepcopy(visited)))
                done = False
                for list1 in walk_lists:
                    if done:
                        break
                    for list2 in walk_lists:
                        if done:
                            break
                        if list1 != list2 and list1[len(list1) - 1] == list2[1] and len(list1) <= 4:
                            self.remove_single_loop(list1)
                            done = True

    def remove_loops(self, input_image, origin):
        """This method removes one side of a loop in the graph (a cycle). The
        side of the loop with the least density will be removed."""
        for node in self.nodes:
            if node.get_num_edges() >= 3: # This is a trinary connection
                walk_lists = list()
                visited = list()
                visited.append(node.get_location())
                for edge in node.get_edges():
                    walk_lists.append(self.walk_until_trinary(edge, deepcopy(visited)))

                done = False
                for list1 in walk_lists:
                    if done:
                        break
                    for list2 in walk_lists:
                        if done:
                            break

                        # This removed a loop along the backbone that ends in a different trinary node.
                        if list1 != list2 and list1[len(list1) - 1] == list2[len(list2) - 1]:
                            density1 = calculate_density(list1, input_image, origin)
                            density2 = calculate_density(list2, input_image, origin)
                            if density1 < density2 and len(list1) <= 4:
                                self.remove_pairs(list1[0], list1[1])
                            elif density2 < density1 and len(list2) <= 4:
                                self.remove_pairs(list2[0], list2[1])
                            done = True

    def remove_pairs(self, location_back, location_front):
        """This is a path-walking method used to remove nodes from the graph

        It is a helper function used by other methods in this file. Each call
        to this method will remove a single node from the graph."""
        node1 = self.get_node(location_back)
        node1.edges.remove(location_front)
        node2 = self.get_node(location_front)
        node2.edges.remove(location_back)
        num_edges = node2.get_num_edges()
        if num_edges == 1:
            self.remove_pairs(location_front, node2.get_edges()[0]) # Should only be one left

    def remove_side_chains(self):
        """This method removes sides chains from the graph

        This method looks for nodes that have three or more edges. It then
        calculates if a single trace connects to another node with three or
        more edges. If this is the case then this trace is likely a side chain
        shortcut that should not exist. It is then removed by this method."""
        for node in self.nodes:
            if node.get_num_edges() >= 3: # This is a trinary connection
                visited = list()
                visited.append(node.get_location())
                min1 = 111 # Smallest
                min1_edge = None
                min2 = 999
                for edge in node.get_edges():
                    value = self.walk_graph(edge, deepcopy(visited), 1)
                    if value < min2 and value < min1:
                        min2 = min1
                        min1 = value
                        min1_edge = edge
                    elif value < min2:
                        min2 = value
                if min1 <= 3 <= min2 - min1: # remove from graph
                    self.remove_pairs(node.get_location(), min1_edge)

    def remove_empty_nodes(self):
        """This function removes empty nodes in the graph. (Garbage collection)"""
        for node in self.nodes:
            if node.get_num_edges() == 0:
                self.nodes.remove(node)

    def print_traces(self, sheet_image, helix_image, offset, pdb_file):
        """This method prints all traces in the graph to a single .PDB file"""
        writer = open(pdb_file, 'w')
        already_written = list()
        set_of_traces = list()
        for node in self.nodes:
            if node.get_num_edges() == 1 or node.get_num_edges() > 2:
                for edge in node.get_edges():
                    trace = list()
                    trace.append(node.get_location())
                    previous = node
                    try:
                        cur = self.get_node(edge)
                        while True:
                            trace.append(cur.get_location())
                            if cur.get_num_edges() != 2:
                                break
                            else:
                                temp_node = cur
                                if cur.get_edges()[0] == previous.get_location():
                                    cur = self.get_node(cur.get_edges()[1])
                                else:
                                    cur = self.get_node(cur.get_edges()[0])
                                previous = temp_node
                    except ValueError as e:
                        print(e)

                    representation = repr(trace[0]) + repr(trace[len(trace) - 1])
                    if representation not in already_written:
                        set_of_traces.append(trace)
                        already_written.append(repr(trace[len(trace) - 1]) + repr(trace[0]))

        helix_traces = list()
        helix_chains = list()
        sheet_traces = list()
        sheet_chains = list()
        cur_chain = 0
        counter = 0
        for trace in set_of_traces:
            cur_chain += 1
            cur_helix = None
            cur_sheet = None
            for ca in trace:
                ca_orig = [int(ca[2] - offset[2]), int(ca[1] - offset[1]), int(ca[0] - offset[0])]
                if helix_image[ca_orig[0]][ca_orig[1]][ca_orig[2]] > 0:
                    if cur_helix is None:
                        cur_helix = list()
                        cur_helix.append(counter)
                    else:
                        cur_helix.append(counter)
                else:
                    if cur_helix is not None:
                        helix_traces.append(cur_helix)
                        helix_chains.append(cur_chain)
                        cur_helix = None
                if sheet_image[ca_orig[0]][ca_orig[1]][ca_orig[2]] > 0:
                    if cur_sheet is None:
                        cur_sheet = list()
                        cur_sheet.append(counter)
                    else:
                        cur_sheet.append(counter)
                else:
                    if cur_sheet is not None:
                        sheet_traces.append(cur_sheet)
                        sheet_chains.append(cur_chain)
                        cur_sheet = None

                PDB_Reader_Writer.write_single_pdb(file=writer, type='ATOM', chain='A', node=np.array([ca[0],ca[1],ca[2]]), seqnum=(counter + cur_chain))
                counter += 1
            PDB_Reader_Writer.write_single_pdb(file=writer, type='TER')

            if cur_helix is not None: # Fence-Posting
                helix_traces.append(cur_helix)
                helix_chains.append(cur_chain)
            if cur_sheet is not None: # Fence-Posting
                sheet_traces.append(cur_sheet)
                sheet_chains.append(cur_chain)

        for index, trace in enumerate(helix_traces):
            chain = helix_chains[index]
            start = str(helix_traces[index][0] + chain - 1)
            end = str(helix_traces[index][len(helix_traces[index]) - 1] + chain - 1)
            PDB_Reader_Writer.write_single_pdb(file=writer, type='HELIX', chain='A', node_from=start, node_to=end)

        for index, trace in enumerate(sheet_traces):
            chain = sheet_chains[index]
            start = str(sheet_traces[index][0] + chain - 1)
            end = str(sheet_traces[index][len(sheet_traces[index]) - 1] + chain - 1)
            PDB_Reader_Writer.write_single_pdb(file=writer, type='SHEET', chain='A', node_from=start, node_to=end)

        writer.close()

    def refine_backbone(self, backbone_image, origin):
        box_size = np.shape(backbone_image)
        new_backbone = np.zeros(box_size)
        steps = 10
        already_written = list()
        for node in self.nodes:
            for edge in node.get_edges():
                node_location = node.get_location()
                node_location = [int(node_location[2] - origin[2]),
                                 int(node_location[1] - origin[1]),
                                 int(node_location[0] - origin[0])] # Reverse it
                edge_location = [int(edge[2] - origin[2]),
                                 int(edge[1] - origin[1]),
                                 int(edge[0] - origin[0])] # Reverse it
                representation = repr(edge_location) + repr(node_location)
                if representation not in already_written:
                    already_written.append(repr(node_location) + repr(edge_location))
                    midpoints = list()
                    x_step = (node_location[0] - edge_location[0]) / steps
                    y_step = (node_location[1] - edge_location[1]) / steps
                    z_step = (node_location[2] - edge_location[2]) / steps
                    for index in range(steps + 1):
                        midpoints.append([edge_location[0] + x_step * index,
                                          edge_location[1] + y_step * index,
                                          edge_location[2] + z_step * index])
                    for z in range(int(edge_location[2]) - 4, int(edge_location[2]) + 4): # 4 is kinda arbitrary here
                        for y in range(int(edge_location[1]) - 4, int(edge_location[1]) + 4):
                            for x in range(int(edge_location[0]) - 4, int(edge_location[0]) + 4):
                                placed = False
                                for index in range(len(midpoints)):
                                    if (box_size[2] > z >= 0 <= y < box_size[1] and 0 <= x < box_size[0] and
                                        distance((x, y, z), midpoints[index]) <= 2
                                            and not placed):
                                        placed = True
                                        new_backbone[x][y][z] = backbone_image[x][y][z]
        new_backbone = np.array(new_backbone, dtype=np.float32)
        return new_backbone

    def print_graph(self, pdb_file):
        """This method prints the graph as a list of edges

        There is no ordering in the graph. It is merely a tool for getting a
        feel for each connection between each Ca atom."""
        writer = open(pdb_file, 'w')
        already_written = list()
        counter = 1
        for node in self.nodes:
            for edge in node.get_edges():
                node_location = node.get_location()
                representation = repr(edge) + repr(node_location)
                if representation not in already_written:
                    PDB_Reader_Writer.write_single_pdb(file=writer, type='ATOM', chain='A', node=np.array([node_location[0],node_location[1],node_location[2]]), seqnum=counter)
                    PDB_Reader_Writer.write_single_pdb(file=writer, type='ATOM', chain='A', node=np.array([edge[0],edge[1],edge[2]]), seqnum=(counter + 1))
                    counter += 3
                    already_written.append(repr(node_location) + repr(edge))
        writer.close()


def calculate_density(walk_list, full_image, origin):
    """Calculates the density of a trace

    This is used to determine which trace should be removed when there are
    more than two traces coming out of a node"""
    box_size = np.shape(full_image)
    steps = 10
    total_density = 0
    number_of_points = 0
    for i in range(len(walk_list) - 1):
        start_point_trans = [walk_list[i][2] - origin[2],
                             walk_list[i][1] - origin[1],
                             walk_list[i][0] - origin[0]]
        end_point_trans = [walk_list[i + 1][2] - origin[2],
                           walk_list[i + 1][1] - origin[1],
                           walk_list[i + 1][0] - origin[0]]
        midpoints = list()
        x_step = (end_point_trans[0] - start_point_trans[0]) / steps
        y_step = (end_point_trans[1] - start_point_trans[1]) / steps
        z_step = (end_point_trans[2] - start_point_trans[2]) / steps
        for j in range(steps + 1):
            midpoints.append([start_point_trans[0] + x_step * j,
                              start_point_trans[1] + y_step * j,
                              start_point_trans[2] + z_step * j])
        for z in range(int(start_point_trans[2]) - 4, int(start_point_trans[2]) + 4): # 4 is kinda arbitrary here
            for y in range(int(start_point_trans[1]) - 4, int(start_point_trans[1]) + 4):
                for x in range(int(start_point_trans[0]) - 4, int(start_point_trans[0]) + 4):
                    placed = False
                    for j in range(len(midpoints)):
                        if (box_size[2] > z >= 0 <= y < box_size[1] and 0 <= x < box_size[0] and
                                    distance((x, y, z), midpoints[j]) <= 1
                            and not placed):
                            placed = True
                            total_density += full_image[x][y][z]
                            number_of_points += 1

    return total_density / number_of_points


def make_graph(pdb_file):
    """This function is not part of the Graph or Node class but it is
    essentially a wrapper method. This method is called to turn an input .PDB
    file into a Graph representation for later processing. This allows the
    prediction step to be separated from the post-processing step. They do not
    have to be run concurrently."""
    pdb_file = open(pdb_file, 'r')
    graph = Graph()
    previous_location = None
    cur_index = -1
    for line in pdb_file:
        if line.startswith("ATOM"):
            index = PDB_Reader_Writer.read_single_pdb_line(type='ATOM INDEX', line=line)
            x, y, z = PDB_Reader_Writer.read_single_pdb_line(type='ATOM', line=line)
            if index == cur_index + 1:
                previous = graph.get_node(previous_location)
                previous.add_edge([x, y, z])
                if graph.contains_location([x, y, z]):
                    node = graph.get_node([x, y, z])
                    node.add_edge(previous_location)
                else:
                    new_node = Node([x, y, z])
                    new_node.add_edge(previous.get_location())
                    graph.add_node(new_node)
            else: # new chain
                if not graph.contains_location([x, y, z]):
                    new_node = Node([x, y, z])
                    graph.add_node(new_node)
            cur_index = index
            previous_location = [x, y, z] # Update for next go-around
    return graph
