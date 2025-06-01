DOI 10.5281/zenodo.15569206

RNA Grammar Variational Autoencoder

The codes in this repository are provided for the generation and optimization of RNA sequences for different optimization critera. 

The required packages are provided in the requirements.txt file. The packages are compatible with Python 3.6.6.


Creating datasets:

python make_RNA_dataset_grammar_py3.py

Training:

python train_RNA.py

Sampling:

python encode_decode_RNA.py

Bayesian optimization:

By going to the folder Theano-master, the modified version of theano could be insalled by typing

python setup.py install

The Bayesian optimization for each of the optimization problems could be done by going to the folder of each of the optimizations and following the provided steps below:

1 - Go to the folder "latent_features_and_targets_grammar", and type:

python generate_latent_features_and_targets.py
2 - Go to the folders "simulation1/grammar/" in the folder of each of the optimizations and type:

nohup python run_bo.py &
Then, repeat the abovementioned command for all the simulation folders (simulation2,...,simulation10).

3 - Go to the folders "simulation1/grammar/results" in the folder of each of the optimizations and type:

python combine_results.py

4 - Going to the main folder of each of the optimization folders and typing:
python combine_results.py
python select_results.py
