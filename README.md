This is the description (documentation) of the simulator.

# 1. Mechanism

We will simulate according to exactly the neural model assumed (the mixed-layered network model).

This simulator is only used for testing the training code (training algorithm, specifically, the stochastic gradient descent used in the mixed-layered neural network model). This simulator is not simulating the true underlying biological mechanism.

I have prepared the tissue hierarchical regulation code (the generation), but now it's not the time for testing that, as we need to test the convergence of the code without this tissue hierarchy first.


# 2. Data generation

The current dir is for code of the simulator, and simulated dataset is in "../simulation_data/" and "../simulation_data_reformat/". The first folder is for the natural output of the simulated data, while the "reformat" folder is for the simulated data in standard format (that the training program recognizes).

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
```

so make sure you have the access to create them when running the program.


# 3. How to use


xxx
