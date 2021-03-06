# 1. Algorithms

## 1.1 Main simulation

We will simulate according to exactly the neural model assumed (the mixed-layered network model). Specifically, there are three pathways from the SNP and batch variables to the gene expression level, the **_cis_- regulation** (each gene has it's own _cis_- regulators), **hidden cell factor**, and the **nonlinear batch effect**. The genotype and batch variables are from uniformly random drawing (0-1), and the five parts of parameters are from Gaussian drawing with Spike and Slab sparsity prior. It's optional that tissue hierarchy can be used as a prior on tissue-specific parameters (_cis_- regulation, cell factor regulation). The gene expression levels (and other hidden variables) are calculated from the simulated variables and the coefficients with the given graphical model.

The source data to be simulated (or generated afterwards) are:

1. genotype
2. batch variables
3. gene expression levels

The model (coefficients) to be simulated is:

1. _cis_- regulation (with hierarchy prior or not)
2. batch to batch hidden regulation
3. batch hidden to gene expression regulation
4. genotype to cell factor regulation
5. cell factor to gene expression regulation (with hierarchy prior or not)

In next step, we can probably build a better simulator with the true genotype but simulated model (the principle from Itsik: mingle the simulation with real data). But I think the current simulation of these variables (genotype and batch variables) won't affect us to evaluate our code of training program.

Other notes of the simulator:

1. We have all the coefficients one more item, which is the intercept of the regression: cis-SNP, SNP-cell factor, cell factor-gene, batch-batch factor, batch factor-gene
2. we might always simulate in the single chromosome mode
3. we might always simulate a full tensor, and then pick up the training set and testing set afterwards
4. the names of individuals/genes/tissues are indices coded in numbers from 0
5. there is no "gene_xymt.txt" information from simulation, so it's an empty file for the training program


## 1.2 Hierarchy regulator

The hierarchy regulation part assumes a given hierarchy, and generate all the nodes from root, with a Gaussian with mean as parent node and variance as the corresponding branch length.

**_Note:_** More discussions will be followed when the code is working in the program.


## 1.3 Partitioner

There is also a partitioner. This is used to partition the full expression tensor (individual x tissue x gene) into traning set and testing set. As the simulation program bases all individual/tissue/gene variables as 0-xxx, here we only specify the range of them, and their indices will be used as their IDs. We partition the samples in each tissue into training set and testing set, according to a given ratio.

## 1.4 More

1. This simulator is only used for testing the training code (training algorithm, specifically, the stochastic gradient descent used in the mixed-layered neural network model). This simulator is not simulating the true underlying biological mechanism (though the model is indeed a good abstract of the understood biology).
2. [**_Jan.12, 2016_**] I have prepared the tissue hierarchical regulation code (the generation), but now it's not the time for testing that, as we need to test the convergence of the code without this tissue hierarchy first.


# 2. Data generation

The current dir is for code of the simulator, and simulated dataset is in **"../simulation_data/"**, **"../simulation_data_reformat/"** and **"../simulation_para/"**. The first folder is for the natural output of the simulated data (that was designed with the first version of this simulator); the "reformat" folder is for the simulated data in standard format (that the training program recognizes; the training program will only load the source data from this folder); and the "para" folder is for the parameters in this model (all the parameters along three pathways), that are in standard format the training program recognizes (these parameters might be used as the model initialization in the training program).

Please see [the training program repo](https://github.com/morrisyoung/eQTL_cplusplus) for the standard data format (source data, model parameters) the training program recognizes.


Note, the simulator will generate the following folders if they are not there:

```
../simulation_data
../simulation_data/SNP_par
../simulation_data/cell_par_cell_gene
../simulation_data/gene

../simulation_data_reformat
../simulation_data_reformat/SNP_par
../simulation_data_reformat/cell_par_cell_gene
../simulation_data_reformat/genotype
../simulation_data_reformat/genotype/chr1

../simulation_para
../simulation_para/para_init_cis_gene
../simulation_para/para_init_cellenv_gene
```
(there are still some redundency in the second folder, as we won't load the model parameters from there, but we still save a copy of the model parameters there)

so make sure you have the access to create them when running the program.


# 3. Dependency


1. Numpy [http://www.numpy.org/]
2. Scipy [http://www.scipy.org/]


# 4. How to use

1. the simulator: python simulator.py (--> ../simulation_data/  ../simulation_data_reformat/  ../simulation_para/)
2. the error_adder: python error_adder.py (--> ../simulation_para_with_error/)
3. the partitioner (training set and testing set; need to set the length of three dimensions): python partition_train_test.py (--> ../simulation_data_reformat/list_samples_train.txt and list_samples_train.txt)



# 5. Logs

Last updated on **_Mar.22, 2016_**.

Last updated on **_Jan.13, 2016_**.




