## Input: SNPs (independent sites), 0/1 binary array
## Output: expression level for some genes
## function: generate the coefficients according to the specified graphical model (the neural model in a Bayesian approach), and generate the corresponding expression level for ALL genes


import numpy as np




"""
n_SNP = 10000
n_eQTL = 10
index_eQTL = []
power_eQTL = []

n_TF = 50
n_Tissue = 13
n_TFSNP = 20  ## or randomly draw from a Gaussian?
power_beforeTF = []
power_afterTF = []
"""




SNP_list = []





"""
## function: randomly generate the genotypes (minor allele frequency) of n individuals for n_SNP independent sites, value in [0, 0.5, 1]
n = 185
n_SNP = 10000


if __name__ == '__main__':

	file = open("genotype.data", "w")


	for i in range(n):

		array = np.random.randint(3, size=n_SNP)
		array = array * 0.5
		for freq in array:
			file.write(str(freq) + ' ')
		file.write('\n')

	file.close()
"""










if __name__ == '__main__':





	"""
	###==============###
	## generate the positions of associated SNPs with the current gene, fill in index_eQTL
	index_eQTL = np.sort( np.random.permutation(n_SNP)[0: n_eQTL] )  ## sorted index
	power_eQTL = 0.5 * np.random.random_sample((n_eQTL,)) + 0.5



	###==============###
	## generate the TF (cell state function) part --> tissue specificity and consistency
	for i in range(n_Tissue):
		array = np.zeros((n_TF, n_SNP))
		power_TF.append(array)
	for i in range(n_TF):
		## generate the index of associated SNPs, and generate their powers, and sign in matrix across different tissues
		index = np.sort( np.random.permutation(n_SNP)[0: n_TFSNP] )  ## sorted index
		power = 0.5 * np.random.random_sample((n_TFSNP, )) + 0.5

		permutation = np.random.permutation(n_Tissue)

		# keep (n_Tissue - 3) amount of tissues as consistent part
		for j in range(n_Tissue - 3):
			tissue = permutation[j]
			## assign all the power of SNPs into this tissue
			for k in range(len(index)):
				power_TF[tissue][index[k]] = power[k]

		# make (3) tissues show tissue specificity; keep the index as the same, but generate the new power
		for j in range(3):
			power = 0.5 * np.random.random_sample((n_TFSNP, )) + 0.5
			tissue = permutation[n_Tissue - 3 + j]
			## assign all the power of SNPs into this tissue
			for k in range(len(index)):
				power_TF[tissue][index[k]] = power[k]


	## generate the TF (cell state function) part --> the power of TF
	power_afterTF = 0.5 * np.random.random_sample((n_TF, )) + 0.5




	###==============###
	## generate the population part
	"""





	###
	

