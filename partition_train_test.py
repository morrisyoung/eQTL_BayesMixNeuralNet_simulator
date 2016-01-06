## this is used to partition the full expression tensor (individual x tissue x gene) into traning set and testing set
## as the simulation program bases all individual/tissue/gene variables as 0-xxx, here we only specify the range of them
## we partition the samples in each tissue into training set and testing set, according to a given ratio


##================
##==== libraries
##================
import numpy as np


##======================
##==== global variables
##======================
# TODO: parameters tunable
n_individual = 185			# as in GTEx.v.4
n_tissue = 17				# as in GTEx.v.4

ratio = 0.75				# the portion of samples used in training


## targets: list_samples_train.txt, list_samples_test.txt
if __name__ == '__main__':

	file1 = open("../simulation_data_reformat/list_samples_train.txt", 'w')
	file2 = open("../simulation_data_reformat/list_samples_test.txt", 'w')

	individual_list = []
	for i in range(n_individual):
		individual_list.append(i)
	individual_list = np.array(individual_list)

	for i in range(n_tissue):
		file1.write(str(i) + '\t')
		file2.write(str(i) + '\t')

		individual_list = np.random.permutation(individual_list)

		max = int(n_individual * 0.75)
		for j in range(max):
			individual = str(individual_list[j])
			tissue = str(i)
			sample = individual + '-' + tissue
			file1.write(sample + '\t')

		for j in range(max, n_individual):
			individual = str(individual_list[j])
			tissue = str(i)
			sample = individual + '-' + tissue
			file2.write(sample + '\t')

		file1.write('\n')
		file2.write('\n')

	file1.close()
	file2.close()


