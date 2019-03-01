from copy import deepcopy
import numpy
import math
import mrcfile
import json


def update_paths(paths):
    paths['normalized_map'] = paths['output'] + 'normalized_map.mrc'


def execute(paths):
    threshold = get_threshold(paths)

    experimental_map = mrcfile.open(paths['cleaned_map'], mode='r')
    experimental_data = deepcopy(experimental_map.data)

    # Remove low valued data and translate the higher values down to zero.
    experimental_data[experimental_data < 0] = 0
    # experimental_data = percentile_filter(experimental_data, numpy.shape(experimental_data), 5)

    # Change all values < threshold to 0
    experimental_data[experimental_data < threshold] = 0
    # translate data to have min = 0
    experimental_data[experimental_data > 0] -= threshold

    # normalize data with percentile value
    percentile = numpy.percentile(experimental_data[numpy.nonzero(experimental_data)], 60)
    experimental_data /= percentile

    # Get rid of the very high-intensity voxels by setting them to 95-percentile
    percentile_98 = numpy.percentile(experimental_data[numpy.nonzero(experimental_data)], 98)
    experimental_data[experimental_data > percentile_98] = percentile_98

    # Print the normalized file to disk.
    with mrcfile.new(paths['normalized_map'], overwrite=True) as mrc:
        mrc.set_data(experimental_data)
        mrc.header.origin = experimental_map.header.origin
        mrc.close()


def get_threshold(paths):
    if 'thresholds_file' in paths:
        emdb_id = paths['input'].split('/')[-2]

        with open(paths['thresholds_file']) as f:
            thresholds = json.load(f)

        if emdb_id in thresholds:
            return thresholds[emdb_id]

    with open(paths['threshold']) as f:
        return float(f.readline())

def distance(z1, z2, y1, y2, x1, x2):
    """Calculates Euclidean distance between two points"""
    z_diff = z1 - z2
    y_diff = y1 - y2
    x_diff = x1 - x2
    sum_squares = math.pow(z_diff, 2) + math.pow(y_diff, 2) + math.pow(x_diff, 2)
    return math.sqrt(sum_squares)


def percentile_filter(full_image, box_size, sphere_radius):
    non_zero_values = list()
    percentile_image = numpy.zeros(box_size)
    filtered_image = numpy.zeros(box_size)
    print('box_size[2]: ' + str(box_size[2]))
    for z in range(box_size[2]):
        print('New z=' + str(z))
        for y in range(box_size[1]):
            for x in range(box_size[0]):
                if full_image[x][y][z] > 0:
                    local_points = list()
                    non_zero_values.append(full_image[x][y][z])
                    # Look at the neighbors
                    for z_n in range(-sphere_radius + z, sphere_radius + z):
                        for y_n in range(-sphere_radius + y, sphere_radius + y):
                            for x_n in range(-sphere_radius + x, sphere_radius + x):
                                if (0 <= z_n < box_size[2] and 0 <= y_n < box_size[1] and 0 <= x_n < box_size[0] and
                                        full_image[x_n][y_n][z_n] > 0 and
                                        distance(z, z_n, y, y_n, x, x_n) <= sphere_radius):
                                    local_points.append(full_image[x_n][y_n][z_n])
                    percentile_90 = numpy.percentile(local_points, 90)
                    percentile_image[x][y][z] = percentile_90
    global_percentile_90 = numpy.percentile(non_zero_values, 90)
    for z in range(box_size[2]):
        for y in range(box_size[1]):
            for x in range(box_size[0]):
                if full_image[x][y][z] > 0:
                    filtered_image[x][y][z] = full_image[x][y][z] * global_percentile_90 / percentile_image[x][y][z]
    return numpy.array(filtered_image, dtype=numpy.float32)
