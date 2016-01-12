This is the description (documentation) of the simulator.

# 1. Algorithms

We will simulate according to exactly the neural model assumed (the mixed-layered network model). Specifically, there are three pathways from the SNP and batch variables to the gene expression level, the _cis_- regulation (each gene has it's own _cis_- regulators), cell factor hidden layer, and the nonlinear batch effect.

The hierarchy regulation part assumes a given hierarchy, and generate all the nodes from root, with a Gaussian with mean as parent node and variance as the corresponding branch length.

One thing to mension is that, this simulator is only used for testing the training code (training algorithm, specifically, the stochastic gradient descent used in the mixed-layered neural network model). This simulator is not simulating the true underlying biological mechanism (though the model is indeed a good abstract of the understood biology).

[**_Jan.12, 2016_**] I have prepared the tissue hierarchical regulation code (the generation), but now it's not the time for testing that, as we need to test the convergence of the code without this tissue hierarchy first.


# 2. Data generation

The current dir is for code of the simulator, and simulated dataset is in "../simulation_data/", "../simulation_data_reformat/" and "../simulation_para/". The first folder is for the natural output of the simulated data (that was designed with the first version of this simulator); the "reformat" folder is for the simulated data in standard format (that the training program recognizes); and the "para" folder is for the parameters in this model (all the parameters along three pathways), that are in standard format the training program recognizes.


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

so make sure you have the access to create them when running the program.


# 3. Dependency


1. Numpy [http://www.numpy.org/]
2. Scipy [http://www.scipy.org/]


# 4. How to use

python simulator.py



# 5. Logs

Last updated on **_Jan.12, 2016_**.


