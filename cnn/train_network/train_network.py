###################################################################################################
# This file is the main script that trains a full Ca/backbone network from start to finish.
# This CNN is actually a series of three CNNs cascaded together. The 1st network does SS
# structure detection. The 2nd network does backbone detection. The 3rd network produces the final
# Ca detection map. The output from each upstream network flows into all downstream networks.
# Moritz, Spencer  November, 2018
###################################################################################################

import tensorflow as tf
import numpy as np
import h5py
import training_parameters as params
import data_manager as dm
import subprocess
from subprocess import Popen
import ca_high_res_network as network

# Create placeholders for the protein data
protein_maps = tf.placeholder('float', name="protein_maps", shape=None)
ss_labels = tf.placeholder('int64', name="ss_labels", shape=None)
backbone_labels = tf.placeholder('int64', name="backbone_labels", shape=None)
ca_labels = tf.placeholder('int64', name="ca_labels", shape=None)

####################################### train_neural_network ######################################
def train_network():
    # SS-Neural Network
    ss_prediction = network.ss_network_definition(protein_maps)
    ss_one_hot_Y = tf.one_hot(ss_labels, 3)
    ss_loss = tf.losses.sigmoid_cross_entropy(multi_class_labels=ss_one_hot_Y, logits=ss_prediction)
    ss_cost = tf.reduce_sum(ss_loss)
    ss_class_predictions = tf.argmax(ss_prediction, 4, name='ss_final_prediction') # 4 indicates the right-most dimention

    # Assign feature map for a-helix and b-sheet
    zeros = tf.zeros(tf.shape(ss_class_predictions))
    ones = tf.ones(tf.shape(ss_class_predictions))
    loops = tf.where(tf.equal(ss_class_predictions, 0), ones, zeros, name='loops_prediction')
    sheet = tf.where(tf.equal(ss_class_predictions, 1), ones, zeros, name='sheet_prediction')
    helix = tf.where(tf.equal(ss_class_predictions, 2), ones, zeros, name='helix_prediction')

    # Backbone-Neural Network
    backbone_prediction = network.backbone_network_definition(protein_maps, sheet, helix)
    backbone_one_hot_Y = tf.one_hot(backbone_labels, 2)
    backbone_loss = tf.losses.sigmoid_cross_entropy(multi_class_labels=backbone_one_hot_Y, logits=backbone_prediction)
    backbone_cost = tf.reduce_sum(backbone_loss)
    backbone_class_predictions = tf.argmax(backbone_prediction, 4, name='backbone_final_prediction') # 4 indicates the right-most dimention

    # CA-Neural Network
    ca_prediction = network.ca_network_definition(protein_maps, sheet, helix, backbone_class_predictions)
    ca_one_hot_Y = tf.one_hot(ca_labels, len(params.classes))
    ca_loss = tf.losses.sigmoid_cross_entropy(multi_class_labels=ca_one_hot_Y, logits=ca_prediction)
    ca_cost = tf.reduce_sum(ca_loss)
    ca_class_predictions = tf.argmax(ca_prediction, 4, name='ca_final_prediction') # 4 indicates the right-most dimention

    # Back-Propogation (loss calculation)
    total_cost = ss_cost + backbone_cost + ca_cost
    optimizer = tf.train.AdamOptimizer(params.learning_rate).minimize(total_cost)

    ### This section is all for Tensorboard metrics collection.
    flattened_labels = tf.reshape(ca_labels, shape=[-1])
    flattened_predictions = tf.reshape(ca_class_predictions, shape=[-1])
    accuracy_run_list, accuracy_opp_list = create_class_accuracy_tensors(ss_labels, ss_class_predictions, backbone_labels, backbone_class_predictions, ca_labels, ca_class_predictions)
    ones = tf.ones(tf.shape(protein_maps))
    zeros = tf.zeros(tf.shape(protein_maps))
    all_class_mask = tf.where(tf.not_equal(protein_maps, 0), ones, zeros)
    flattened_all_class_mask = tf.reshape(all_class_mask, shape=[-1])
    class_accuracy_mean, class_accuracy_mean_op = tf.metrics.mean_per_class_accuracy(labels=flattened_labels, predictions=flattened_predictions, num_classes=len(params.classes), weights=flattened_all_class_mask)
    mean_squared_error, mean_squared_error_op = tf.metrics.mean_squared_error(labels=flattened_labels, predictions=flattened_predictions, weights=flattened_all_class_mask)
    voxelwise_accuracy, voxelwise_accuracy_op = tf.metrics.accuracy(labels=flattened_labels, predictions=flattened_predictions, weights=flattened_all_class_mask)
    class_accuracy_mean_summary = tf.summary.scalar('Train Mean Class Accuracy', class_accuracy_mean)
    test_class_accuracy_mean_summary = tf.summary.scalar('Test Mean Class Accuracy', class_accuracy_mean)
    mean_squared_error_summary = tf.summary.scalar('Train Mean Squared Error', mean_squared_error)
    train_voxelwise_accuracy_summary = tf.summary.scalar('Train Voxelwise Accuracy', voxelwise_accuracy)
    test_voxelwise_accuracy_summary = tf.summary.scalar('Test Voxelwise Accuracy', voxelwise_accuracy)
    accuracy_opp_list.append(class_accuracy_mean_op)
    accuracy_run_list.append(test_class_accuracy_mean_summary)
    accuracy_opp_list.append(voxelwise_accuracy_op)
    accuracy_run_list.append(test_voxelwise_accuracy_summary)
    ### End of Tensorboard metrics collection.

    with tf.Session() as sess:

        tensorboard_writer = tf.summary.FileWriter(params.tensorboard_path, sess.graph) # TensorBoard Write objects

        saver = tf.train.Saver(max_to_keep=10) # Only save a single model at a time. (Overwrite old one)
        if params.restore_model: # Possibly restore a saved model.
            saver.restore(sess, params.saved_module_path + '/saved_model.ckpt')
        else: # Otherwise start from scrath
            sess.run(tf.initialize_all_variables())

        for epoch in range(params.start_epoc, params.end_epochs):
            epoch_loss = 0
            sess.run(tf.local_variables_initializer()) # Here for the tf.metrics

            # For each epoch, feed all train data to the network in a series of batches
            total_batches = trainManager.total_batches()
            for batch in range(total_batches):
                print('Batch ' + str(batch + 1) + '/' + str(total_batches))
                Protein_maps, SS_Labels, Backbone_Labels, CA_Labels = trainManager.next_batch()
                _, c, _, _, _ = sess.run([optimizer, total_cost, class_accuracy_mean_op, mean_squared_error_op, voxelwise_accuracy_op], feed_dict={protein_maps:Protein_maps, ss_labels:SS_Labels, backbone_labels:Backbone_Labels, ca_labels:CA_Labels})
                epoch_loss += c

            # Write the cost and mean_class_accuracy to tensorboard
            training_class_accuracy_mean, training_mean_squared_error, train_voxelwise_accuracy = sess.run([class_accuracy_mean_summary, mean_squared_error_summary, train_voxelwise_accuracy_summary])
            tensorboard_writer.add_summary(training_class_accuracy_mean, epoch)
            tensorboard_writer.add_summary(training_mean_squared_error, epoch)
            tensorboard_writer.add_summary(train_voxelwise_accuracy, epoch)
            
            print('Epoch', epoch, 'completed out of', params.end_epochs,'loss:',epoch_loss)

            sess.run(tf.local_variables_initializer()) # Reset local variables for test run
            # For each epoch, feed all train data to the network.
            total_batches = testManager.total_batches()
            for batch in range(testManager.total_batches()):
                print('Batch ' + str(batch + 1) + '/' + str(total_batches))
                Protein_maps, SS_Labels, Backbone_Labels, CA_Labels = testManager.next_batch()
                results = sess.run(accuracy_opp_list, feed_dict={protein_maps:Protein_maps, ss_labels:SS_Labels, backbone_labels:Backbone_Labels, ca_labels:CA_Labels})
            evaluate_category_accuracy(session=sess, writer=tensorboard_writer, run_list=accuracy_run_list, epoch=epoch)
            save_path = saver.save(sess, params.saved_module_path + '/saved_model.ckpt') # Save the model after each epoch

def create_class_accuracy_tensors(ss_labels, ss_predictions, backbone_labels, backbone_predictions, ca_labels, ca_predictions):
    accuracy_opp_list = []
    accuracy_run_list = []
    ones = tf.ones(tf.shape(ss_labels))
    zeros = tf.zeros(tf.shape(ss_labels))

    # CA-Metrics
    ca_flattened_labels = tf.reshape(ca_labels, shape=[-1])
    ca_flattened_predictions = tf.reshape(ca_predictions, shape=[-1])
    for category in range(len(params.classes)):
        ca_flattened_category_mask = tf.reshape(tf.where(tf.equal(ca_labels, category), ones, zeros), shape=[-1])
        accuracy_run, accuracy_opp = tf.metrics.accuracy(labels=ca_flattened_labels, predictions=ca_flattened_predictions, weights=ca_flattened_category_mask)
        accuracy_opp_list.append(accuracy_opp)
        accuracy_run_list.append(tf.summary.scalar('Test ' + params.classes[category] + ' Class Accuracy', accuracy_run))

    # SS-Metrics
    ss_flattened_labels = tf.reshape(ss_labels, shape=[-1])
    ss_flattened_predictions = tf.reshape(ss_predictions, shape=[-1])
    helix_cats = ['Loops', 'Sheets', 'Helix ']
    for index in range(3):
        ss_flattened_category_mask = tf.reshape(tf.where(tf.equal(ss_labels, index), ones, zeros), shape=[-1])
        accuracy_run, accuracy_opp = tf.metrics.accuracy(labels=ss_flattened_labels, predictions=ss_flattened_predictions, weights=ss_flattened_category_mask)
        accuracy_opp_list.append(accuracy_opp)
        accuracy_run_list.append(tf.summary.scalar('Test ' + helix_cats[index] + ' Class Accuracy', accuracy_run))

    # Backbone-Metrics
    backbone_flattened_labels = tf.reshape(backbone_labels, shape=[-1])
    backbone_flattened_predictions = tf.reshape(backbone_predictions, shape=[-1])
    backbone_cats = ['Sidechain', 'Backbone']
    for index in range(2):
        backbone_flattened_category_mask = tf.reshape(tf.where(tf.equal(backbone_labels, index), ones, zeros), shape=[-1])
        accuracy_run, accuracy_opp = tf.metrics.accuracy(labels=backbone_flattened_labels, predictions=backbone_flattened_predictions, weights=backbone_flattened_category_mask)
        accuracy_opp_list.append(accuracy_opp)
        accuracy_run_list.append(tf.summary.scalar('Test ' + backbone_cats[index] + ' Class Accuracy', accuracy_run))

    return accuracy_run_list, accuracy_opp_list

def evaluate_category_accuracy(session, writer, run_list, epoch):
    results = session.run(run_list)
    for category in range(len(params.classes) + 2 + 2 + 2): # The extra two are for the mean class accuracy and voxelwise accuracy
        writer.add_summary(results[category], epoch)

# Concatinates the protein weight masks to produce an extra dimention suitable for the one-hot-encoding output
def concatinate_masks(weight_mask):
    weight_mask_extended = tf.expand_dims(weight_mask, 4)
    list_of_masks = []
    for index in range(len(params.classes)):
        list_of_masks.append(weight_mask_extended)
    return tf.concat(list_of_masks, 4)

############################################ main #################################################
trainManager = dm.DataManager()
trainManager.openHDF5(params.data_path + params.model + '_train_data.hdf5')
trainManager.set_batch_size(batch_size=params.batch_size)
testManager = dm.DataManager()
testManager.openHDF5(params.data_path + params.model + '_test_data.hdf5')
testManager.set_batch_size(batch_size=params.batch_size)

# This launches Tensorboard for monitoring the traning.
subprocess.Popen('tensorboard --logdir=' + params.tensorboard_path)

###### Train Model ######
train_network()
#########################

# Cleanup
del trainManager
del testManager