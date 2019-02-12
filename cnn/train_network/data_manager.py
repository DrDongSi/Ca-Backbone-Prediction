###################################################################################################
# This file manages the input data for the neural network. This class opens an HDF5 file that
# contains the properly organized train and test data. It then dispenses that data to the neural
# network in batch sizes for training and testing.
# Moritz, Spencer  November, 2018
###################################################################################################

import numpy as np
import h5py
import math

class DataManager:
    def __init__(self):
        self.file_handle = None
        self.protein_maps = []
        self.ss_labels = []
        self.backbone_labels = []
        self.ca_labels = []
        self.protein_ids = []
        self.number_proteins = 0
        self.active_protein = 0
        self.boxSize = 0
        self.batch_size = 0

    def __del__(self):
        self.file_handle.close()

    # Open the HDF5 file from the saved location on disk. Save handels to each piece of data.
    def openHDF5(self, filename):
        self.file_handle = h5py.File(filename, 'r')
        self.protein_maps = self.file_handle['protein_maps']
        self.ss_labels = self.file_handle['ss_labels']
        self.backbone_labels = self.file_handle['backbone_labels']
        self.ca_labels = self.file_handle['ca_labels']
        self.protein_ids  = self.file_handle['protein_ids']
        self.number_proteins = self.protein_maps.shape[0]
        self.boxSize = self.protein_maps.shape[1]

    def set_batch_size(self, batch_size):
        self.batch_size = batch_size

    def total_batches(self):
        return int(math.ceil(self.number_proteins / self.batch_size))

    # Grabs the next batch of data and returns it to caller for training in the CNN.
    def next_batch(self):
        num_proteins = self.batch_size if self.batch_size <= self.number_proteins - self.active_protein else self.number_proteins - self.active_protein
        protein_maps = np.ndarray((num_proteins, self.boxSize, self.boxSize, self.boxSize))
        ss_labels = np.ndarray((num_proteins, self.boxSize, self.boxSize, self.boxSize))
        backbone_labels = np.ndarray((num_proteins, self.boxSize, self.boxSize, self.boxSize))
        ca_labels = np.ndarray((num_proteins, self.boxSize, self.boxSize, self.boxSize))
        for index in range(num_proteins):
            protein_maps[index]  = self.protein_maps[self.active_protein][:][:][:]
            ss_labels[index] = self.ss_labels[self.active_protein][:][:][:]
            backbone_labels[index] = self.backbone_labels[self.active_protein][:][:][:]
            ca_labels[index] = self.ca_labels[self.active_protein][:][:][:]
            self.active_protein = (self.active_protein + 1) % self.number_proteins
        return protein_maps, ss_labels, backbone_labels, ca_labels

