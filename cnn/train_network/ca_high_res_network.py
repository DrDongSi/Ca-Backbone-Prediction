###################################################################################################
# This file contains the definition of the Ca/backbone prediction cascaded neural network. The
# three networks are divided into different functions.
# Moritz, Spencer  November, 2018
###################################################################################################

import tensorflow as tf
import training_parameters as params

########################################## Leaky Relu #############################################
def lrelu(x,alpha=0.1):
    return tf.maximum(alpha*x,x)

############################### Ca atom prediction CNN definition #################################
def ca_network_definition(protein_maps, sheet, helix, backbone_predictions):
    protein_maps = tf.reshape(protein_maps, shape=[-1, params.box_size, params.box_size, params.box_size, 1])
    sheet = tf.reshape(sheet, shape=[-1, params.box_size, params.box_size, params.box_size, 1])
    helix = tf.reshape(helix, shape=[-1, params.box_size, params.box_size, params.box_size, 1])
    backbone_predictions = tf.reshape(backbone_predictions, shape=[-1, params.box_size, params.box_size, params.box_size, 1])
    x = tf.concat([protein_maps, tf.cast(sheet, dtype=tf.float32), tf.cast(helix, dtype=tf.float32), tf.cast(backbone_predictions, dtype=tf.float32)], 4)
    conv1 = tf.layers.conv3d(x, filters=32, kernel_size=(4,4,4), strides=(1,1,1), name='ca_conv1', padding='SAME', use_bias=True)
    conv2 = tf.layers.conv3d(conv1, filters=64, kernel_size=(4,4,4), strides=(1,1,1), padding='SAME', dilation_rate=(2,2,2), use_bias=True, activation=lrelu, name='ca_conv2')
    conv3 = tf.layers.conv3d(conv2, filters=64, kernel_size=(4,4,4), strides=(1,1,1), padding='SAME', dilation_rate=(2,2,2), use_bias=True, activation=lrelu, name='ca_conv3')
    return tf.layers.conv3d(conv3, filters=2, kernel_size=(4,4,4), strides=(1,1,1), name='ca_logits', padding='SAME', use_bias=True)

######################### backbone structure prediction CNN definition ############################
def backbone_network_definition(protein_maps, sheet, helix):
    protein_maps = tf.reshape(protein_maps, shape=[-1, params.box_size, params.box_size, params.box_size, 1])
    sheet = tf.reshape(sheet, shape=[-1, params.box_size, params.box_size, params.box_size, 1])
    helix = tf.reshape(helix, shape=[-1, params.box_size, params.box_size, params.box_size, 1])
    x = tf.concat([protein_maps, tf.cast(sheet, dtype=tf.float32), tf.cast(helix, dtype=tf.float32)], 4)
    conv1 = tf.layers.conv3d(x, filters=32, kernel_size=(4,4,4), strides=(1,1,1), name='backbone_conv1', padding='SAME', use_bias=True)
    conv2 = tf.layers.conv3d(conv1, filters=64, kernel_size=(4,4,4), strides=(1,1,1), padding='SAME', dilation_rate=(2,2,2), use_bias=True, activation=lrelu, name='backbone_conv2')
    conv3 = tf.layers.conv3d(conv2, filters=64, kernel_size=(4,4,4), strides=(1,1,1), padding='SAME', dilation_rate=(2,2,2), use_bias=True, activation=lrelu, name='backbone_conv3')
    return tf.layers.conv3d(conv3, filters=2, kernel_size=(4,4,4), strides=(1,1,1), name='backbone_logits', padding='SAME', use_bias=True)

######################### secondary structure prediction CNN definition ###########################
def ss_network_definition(x):
    x = tf.reshape(x, shape=[-1, params.box_size, params.box_size, params.box_size, 1])
    conv1 = tf.layers.conv3d(x, filters=32, kernel_size=(5,5,5), strides=(1,1,1), name='ss_conv1', padding='SAME', use_bias=True)
    conv2 = tf.layers.conv3d(conv1, filters=64, kernel_size=(4,4,4), strides=(1,1,1), padding='SAME', dilation_rate=(2,2,2), use_bias=True, activation=lrelu, name='ss_conv2')
    conv3 = tf.layers.conv3d(conv2, filters=64, kernel_size=(4,4,4), strides=(1,1,1), padding='SAME', dilation_rate=(2,2,2), use_bias=True, activation=lrelu, name='ss_conv3')
    return tf.layers.conv3d(conv3, filters=3, kernel_size=(5,5,5), strides=(1,1,1), name='ss_logits', padding='SAME', use_bias=True)