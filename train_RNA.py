from __future__ import print_function

import argparse
import os
import h5py
import numpy as np

from models.model_RNA import RNAVAE
#from keras.callbacks import ModelCheckpoint, ReduceLROnPlateau
import tensorflow as tf

import h5py
import RNA_grammar as G

MAX_LEN = 199#137 #64 # 277
'''
LATENT = 25
EPOCHS = 50
BATCH = 600
'''
rules = G.gram.split('\n')
DIM = len(rules)
LATENT = 10 # 2
EPOCHS = 100
BATCH = 100#500

NCHARS = len(G.GCFG.productions())



print("rules")
print(DIM)

def get_arguments():
    parser = argparse.ArgumentParser(description='Molecular autoencoder network')
    parser.add_argument('--epochs', type=int, metavar='N', default=EPOCHS,
                        help='Number of epochs to run during training.')
    parser.add_argument('--latent_dim', type=int, metavar='N', default=LATENT,
                        help='Dimensionality of the latent representation.')
    return parser.parse_args()

def main():

    # 0. load dataset
    h5f = h5py.File('data/RNA_grammar_SLF_dataset_uniform_prob_0_100k_Py3.h5', 'r')
    data = h5f['data'][:]
    h5f.close()
    
    
    
    # 1. split into train/test, we use test set to check reconstruction error and the % of
    # samples from prior p(z) that are valid
    
    print("shape")
    print(np.shape(data))
    XTE = data[0:1000]    #1842
    XTR = data[1000:]
    
    print("type")
    print(type(data))
    print(data[1])
    #data = int(data)
    data2 = data.astype(int)
    print("data2")
    print(data2[1])
    
    print("rules")
    print(DIM)

    # 1. get any arguments and define save file, then create the VAE model
    args = get_arguments()
    #params = {'hidden': 501, 'dense': 435, 'conv1': 9, 'conv2': 9, 'conv3': 10}
    params = {'hidden': 100, 'dense': 100, 'conv1': 9, 'conv2': 9, 'conv3': 10}
    model_save = 'RNA_vae_grammar_SLF_h_100k' + str(params['hidden']) + '_c234_L' + str(args.latent_dim) + '_E' + str(args.epochs) + '_batchB.hdf5'
    model = RNAVAE()

    # 2. if this results file exists already load it
    if os.path.isfile(model_save):
        #model.load(rules, model_save, latent_rep_size = args.latent_dim, hypers = params)
        model.load(NCHARS, model_save, latent_rep_size = args.latent_dim, hypers = params)
    else:
        #model.create(rules, max_length=MAX_LEN, latent_rep_size = args.latent_dim, hypers = params)
        model.create(NCHARS, max_length=MAX_LEN, latent_rep_size = args.latent_dim, hypers = params)

    # 3. only save best model found on a 10% validation set
    '''
    checkpointer = ModelCheckpoint(filepath = model_save,
                                   verbose = 1,
                                   save_best_only = True)

    reduce_lr = ReduceLROnPlateau(monitor = 'val_loss',
                                  factor = 0.2,
                                  patience = 1,
                                  min_lr = 0.0001)
    '''
    checkpointer = tf.keras.callbacks.ModelCheckpoint(filepath = model_save,
                                   verbose = 1,
                                   save_best_only = True)

    reduce_lr = tf.keras.callbacks.ReduceLROnPlateau(monitor = 'val_loss',
                                  factor = 0.2,
                                  patience = 3,
                                  min_lr = 0.0001)

    # 4. fit the vae
    model.autoencoder.fit(
        XTR,
        XTR,
        shuffle = True,
        nb_epoch = args.epochs,
        batch_size = BATCH,
        callbacks = [checkpointer, reduce_lr],
        validation_split = 0.1
    )

if __name__ == '__main__':
    main()
