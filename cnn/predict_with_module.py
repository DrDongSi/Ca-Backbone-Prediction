"""This script runs the Ca/Backbone prediction CNN with a single protein's pre-processed mrc file

The script accomplished the following:
1. Load the pre-processed protein from an .MRC file
2. Restore the saved prediction model from file.
3. Runs the protein through the network using a splitting function for proteins larger than 64^3 in size.
4. Prints the SS prediction maps, backbone prediction map, and Ca prediction map to file.
5. Uses graph theory to improve the final backbone trace.
6. Print the final traces to file.
"""

import tensorflow as tf
import numpy as np
import mrcfile
from copy import deepcopy
import math
import cnn.map_splitter as ms
from collections import deque
import os
import prediction as pre

__author__ = 'Moritz Spencer'


def update_paths(paths):
    paths['loops_confidence'] = paths['output'] + 'loops_confidence.mrc'
    paths['sheet_confidence'] = paths['output'] + 'sheet_confidence.mrc'
    paths['helix_confidence'] = paths['output'] + 'helix_confidence.mrc'
    paths['backbone_confidence'] = paths['output'] + 'backbone_confidence.mrc'
    paths['ca_confidence'] = paths['output'] + 'ca_confidence.mrc'


def execute(paths):
    with pre.prediction.semaphore:
        with tf.Session() as sess:
            module_path = os.path.dirname(os.path.abspath(__file__)) + '/saved_module/5-7A_Full_SS_Combo/'
            saver = tf.train.import_meta_graph(module_path + 'saved_model.ckpt.meta')
            saver.restore(sess, module_path + 'saved_model.ckpt')  # Load the saved CNN.
            normalized_map = mrcfile.open(paths['normalized_map'], mode='r')
            full_image = deepcopy(normalized_map.data)

            manifest = ms.create_manifest(full_image) # Create a 'manifest' to run through the CNN.

            # Tensors to be restored in the CNN. These tensors will hold the final output from each stage.
            graph = tf.get_default_graph()
            x = graph.get_tensor_by_name("protein_maps:0")
            y = graph.get_tensor_by_name("ss_labels:0")
            loops_op = graph.get_tensor_by_name('loops_prediction:0')
            sheet_op = graph.get_tensor_by_name('sheet_prediction:0')
            helix_op = graph.get_tensor_by_name('helix_prediction:0')
            ss_op = graph.get_tensor_by_name('ss_logits/BiasAdd:0')
            backbone_op = graph.get_tensor_by_name('backbone_logits/BiasAdd:0')
            ca_op = graph.get_tensor_by_name('ca_logits/BiasAdd:0')

            # Placeholders for prediction maps that will be output from the CNN.
            loops_image = np.zeros((np.shape(manifest)))
            sheet_image = np.zeros((np.shape(manifest)))
            helix_image = np.zeros((np.shape(manifest)))
            loops_confidence = np.zeros((np.shape(manifest)))
            sheet_confidence = np.zeros((np.shape(manifest)))
            helix_confidence = np.zeros((np.shape(manifest)))
            backbone_image = np.zeros((np.shape(manifest)))
            ca_image = np.zeros((np.shape(manifest)))

            # Run the protein through the CNN and save the output in a local placeholder.
            for index in range(math.ceil(len(manifest) / 10)):
                loops_output, sheet_output, helix_output, ss_output, backbone_output, ca_output = \
                    sess.run([loops_op, sheet_op, helix_op, ss_op, backbone_op, ca_op],
                             feed_dict={
                                 x: manifest[index * 10: (index + 1) * 10],
                                 y: manifest[index * 10: (index + 1) * 10]
                             })
                loops_image[index * 10: (index + 1) * 10] = loops_output
                sheet_image[index * 10: (index + 1) * 10] = sheet_output
                helix_image[index * 10: (index + 1) * 10] = helix_output
                loops_confidence[index * 10: (index + 1) * 10] = ss_output[:, :, :, :, 0]
                sheet_confidence[index * 10: (index + 1) * 10] = ss_output[:, :, :, :, 1]
                helix_confidence[index * 10: (index + 1) * 10] = ss_output[:, :, :, :, 2]
                backbone_image[index * 10: (index + 1) * 10] = np.subtract(backbone_output[:, :, :, :, 1],
                                                                           backbone_output[:, :, :, :, 0])
                ca_image[index * 10: (index + 1) * 10] = np.subtract(ca_output[:, :, :, :, 1], ca_output[:, :, :, :, 0])

            # Add an arbitrary constant for improved viewing in Chimera.
            backbone_image += 4  # Used 4 for sim maps.
            ca_image += 10

            # Reconstruct each 64^3 image into the full protein shape.
            loops_confidence = ms.reconstruct_map(loops_confidence, np.shape(full_image))
            sheet_confidence = ms.reconstruct_map(sheet_confidence, np.shape(full_image))
            helix_confidence = ms.reconstruct_map(helix_confidence, np.shape(full_image))
            backbone_image = ms.reconstruct_map(backbone_image, np.shape(full_image))
            ca_image = ms.reconstruct_map(ca_image, np.shape(full_image))

            # Clean up predicted images by zeroing out space outside in input map.
            input_mask = np.where(full_image > 0, 1, 0)
            loops_confidence = np.where(input_mask == 1, loops_confidence, 0)
            sheet_confidence = np.where(input_mask == 1, sheet_confidence, 0)
            helix_confidence = np.where(input_mask == 1, helix_confidence, 0)
            ss_image = np.stack((loops_confidence, sheet_confidence, helix_confidence), axis=3)
            backbone_image = np.where(input_mask == 1, backbone_image, 0)
            backbone_image[backbone_image < 0] = 0
            ss_image = ss_nearest_neighbor(ss_image, input_mask)  # Post-Processing Step to clean up SS predictions.

            loops_image = np.where(ss_image == 0, 1, 0)
            sheet_image = np.where(ss_image == 1, 1, 0)
            helix_image = np.where(ss_image == 2, 1, 0)
            loops_image = np.where(input_mask == 1, loops_image, 0)
            sheet_image = np.where(input_mask == 1, sheet_image, 0)
            helix_image = np.where(input_mask == 1, helix_image, 0)

            remove_small_chunks(backbone_image)
            ca_image = np.where(input_mask == 1, ca_image, 0)
            ca_image = np.array(ca_image, dtype=np.float32)

            # Print the loops image
            with mrcfile.new(paths['loops_confidence'], overwrite=True) as mrc:
                mrc.set_data(np.array(loops_image, dtype=np.float32))
                mrc.header.origin = normalized_map.header.origin.item(0)
                mrc.update_header_stats()
                mrc.close()

            # Print the sheet image
            with mrcfile.new(paths['sheet_confidence'], overwrite=True) as mrc:
                mrc.set_data(np.array(sheet_image, dtype=np.float32))
                mrc.header.origin = normalized_map.header.origin.item(0)
                mrc.update_header_stats()
                mrc.close()

            # Print the helix image
            with mrcfile.new(paths['helix_confidence'], overwrite=True) as mrc:
                mrc.set_data(np.array(helix_image, dtype=np.float32))
                mrc.header.origin = normalized_map.header.origin.item(0)
                mrc.update_header_stats()
                mrc.close()

            # Print the backbone confidence image
            with mrcfile.new(paths['backbone_confidence'], overwrite=True) as mrc:
                mrc.set_data(backbone_image)
                mrc.header.origin = normalized_map.header.origin.item(0)
                mrc.update_header_stats()
                mrc.close()

            # Print the ca-confidence image
            with mrcfile.new(paths['ca_confidence'], overwrite=True) as mrc:
                mrc.set_data(ca_image)
                mrc.header.origin = normalized_map.header.origin.item(0)
                mrc.update_header_stats()
                mrc.close()

            normalized_map.close()


# Post-Processing step used to remove classification outliers in the secondary structure
# prediction image. This function examines each voxel in the image and reassigns it to
# the secondary structure represented by the weighted average of its neighbors.
def ss_nearest_neighbor(ss_confidence, input_mask):
    box_size = np.shape(input_mask)
    output_prediction = np.zeros(box_size)
    sphere_radius = 2
    for x in range(1, box_size[0] - 1):
        for y in range(1, box_size[1] - 1):
            for z in range(1, box_size[2] - 1):
                if input_mask[x][y][z] > 0:
                    loops_weight = 0
                    sheet_weight = 0
                    helix_weight = 0
                    for z_n in range(-sphere_radius + z, sphere_radius + z):
                        for y_n in range(-sphere_radius + y, sphere_radius + y):
                            for x_n in range(-sphere_radius + x, sphere_radius + x):
                                if (0 <= z_n < box_size[2] and 0 <= y_n < box_size[1] and
                                                0 <= x_n < box_size[0] and input_mask[x_n][y_n][z_n] > 0 and
                                            distance(z, z_n, y, y_n, x, x_n) <= sphere_radius):
                                    loops_weight += ss_confidence[x_n][y_n][z_n][0]
                                    sheet_weight += ss_confidence[x_n][y_n][z_n][1]
                                    helix_weight += ss_confidence[x_n][y_n][z_n][2]
                    if loops_weight > sheet_weight and loops_weight > helix_weight:
                        output_prediction[x][y][z] = 0
                    elif sheet_weight > helix_weight:
                        output_prediction[x][y][z] = 1
                    else:
                        output_prediction[x][y][z] = 2
    return output_prediction


def remove_small_chunks(input_image):
    """A method used to remove small disjoint regions in a 3D image

    This was primarily used to clean up the backbone prediction .MRC file
    however it may not be necessary
    """
    min_chuck_size = 25
    box_size = np.shape(input_image)
    visited = np.zeros(box_size)
    for x in range(1, box_size[0] - 1):
        for y in range(1, box_size[1] - 1):
            for z in range(1, box_size[2] - 1):
                if input_image[x, y, z] > 0 and visited[x, y, z] == 0:
                    chunk_list = list()
                    queue = deque()
                    queue.append([x, y, z])
                    chunk_list.append([x, y, z])
                    while len(queue) > 0:
                        cur_position = queue.popleft()
                        visited[cur_position[0], cur_position[1], cur_position[2]] = 1
                        offsets = [[0, 0, -1], [0, 0, 1], [0, -1, 0], [0, 1, 0], [-1, 0, 0], [1, 0, 0]]
                        for index in range(len(offsets)):
                            x_new = cur_position[0] + offsets[index][0]
                            y_new = cur_position[1] + offsets[index][1]
                            z_new = cur_position[2] + offsets[index][2]
                            if 0 <= x_new < box_size[0] and 0 <= y_new < box_size[1] and 0 <= z_new < box_size[2]:
                                if input_image[x_new, y_new, z_new] > 0 and visited[x_new, y_new, z_new] == 0:
                                    queue.append([x_new, y_new, z_new])
                                    visited[x_new, y_new, z_new] = 1
                                    chunk_list.append([x_new, y_new, z_new])
                    if len(chunk_list) < min_chuck_size:
                        for voxel in chunk_list:
                            input_image[voxel[0], voxel[1], voxel[2]] = 0


def distance(z1, z2, y1, y2, x1, x2):
    """Calculates Euclidean distance between two points"""
    z_diff = z1 - z2
    y_diff = y1 - y2
    x_diff = x1 - x2
    sum_squares = math.pow(z_diff, 2) + math.pow(y_diff, 2) + math.pow(x_diff, 2)
    return math.sqrt(sum_squares)
