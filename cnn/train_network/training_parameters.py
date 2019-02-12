###################################################################################################
# This file contains training parameters used to train the Ca/backbone detection CNN. Changes here
# will affect the traning of the network.
# Moritz, Spencer  November, 2018
###################################################################################################

from enum import Enum
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' # Quiet down the CMD output.

global restore_model
global batch_size
global learning_rate
global number_of_classes
global box_size
global folder_path
global saved_module_path
global tensorboard_path
global start_epoc
global end_epochs
global classes
global data_path
global Topic
global topic
global model

restore_model = False
batch_size = 5
learning_rate = 1e-4
box_size = 64 #64^3

start_epoc = 0
end_epochs = 100

base_data_folder_path = 'D:/data/packaged_data/'
classes = ['Other', 'Ca']
folder_path = './ca_high_res/'
data_path = os.path.join(base_data_folder_path, 'ca_high_res_data/')

model = 'hello_world' # This must match the name of your test/train data.

saved_module_path = os.path.join(folder_path, 'saved_module', model)
tensorboard_path = os.path.join(folder_path, 'tensorboard', model)
