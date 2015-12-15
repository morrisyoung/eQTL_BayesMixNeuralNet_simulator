## Input: SNPs (independent sites), 0/1 binary array
## Output: expression level for some genes
## function: generate the coefficients according to the specified graphical model (the neural model in a Bayesian approach), and generate the corresponding expression level for ALL genes
## notes:
##	1. we should simulate tissue specificity, as the modeling takes consideration of that;
##	2. when we say a factor, it's latent factor condensed from some original variables
##	3. xxx


import numpy as np


##==== global variables
## notes: TODO
##	1. shall we simulate Spike and Slab?
##	2. xxx
##=====================
#==== chromosome
L = 0								# unit: basepair

#==== individual
n_individual = 0

#==== tissue
n_tissue = 0

#==== SNP
n_SNP = 0
SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
SNP_pos_list = []
SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess

#==== gene
n_gene = 0
gene_rep = {individual:{tissue:[], ...}, ...}			# gene expression list for tissue samples for individuals
gene_pos_list = []

#==== SNP gene pos map
pos_map = {}							# {gene:[snp1, snp2], ...}

#==== factor_cell
n_factor_cell = 0
factor_cell_beta_rep = {cell factor:[], ...}			# coefficient lists for cell factors
beta_factor_cell_rep = {tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues

#==== factor_batch (analogous to cell factor pathway)
n_batch = 0							# n_batch = n_batch_individual + n_batch_sample
n_batch_individual = 0
n_batch_sample = 0
#n_factor_batch							# the number of latent batch factors
batch_individual_rep = {individual:[], ...}			# batch variable lists for individuals
batch_tissue_sample_rep = {individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
factor_batch_beta_rep = {batch factor:[], ...}			# coefficient lists for batch factors
beta_factor_batch_rep = {gene:[], ...}				# coefficient lists of batch factors for genes





##==== simulating variables
## notes: TODO
##	1. simulate variables observed and latent excluding those along the generation process
##	2. xxx
##=====================
def simu_genotype():	# variables
	# to fill in: 
	##SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
	##SNP_pos_list = []

"""
## function: randomly generate the genotypes (minor allele frequency) of n individuals for n_SNP independent sites, value in [0, 0.5, 1]
n = 185
n_SNP = 10000

	file = open("genotype.data", "w")


	for i in range(n):

		array = np.random.randint(3, size=n_SNP)
		array = array * 0.5
		for freq in array:
			file.write(str(freq) + ' ')
		file.write('\n')

	file.close()
"""

	return

def simu_genotype_beta():	# beta (cis-); one more coefficient as the constant item in linear regression
	# to fill in:
	##SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	# later on, after "SNP gene pos map"



	return


def simu_gene():
	# to fill in:
	##gene_pos_list = []


	return



def simu_cell_factor():
	# to fill in:
	##factor_cell_beta_rep = {cell factor:[], ...}			# coefficient lists for cell factors
	##beta_factor_cell_rep = {tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues


	return



## simulate batch variables for all individuals, and for all samples in different individuals; then concatenate them when using them
def simu_batch_factor():	# variables and beta
	# to fill in:
	##batch_individual_rep = {individual:[], ...}			# batch variable lists for individuals
	##batch_tissue_sample_rep = {individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
	##factor_batch_beta_rep = {batch factor:[], ...}			# coefficient lists for batch factors
	##beta_factor_batch_rep = {gene:[], ...}


	return











##==== other utilities
## notes: TODO
##	1. xxx
##=====================
def SNP_gene_map():
	# to fill in:
	##pos_map = {}							# {gene:[snp1, snp2], ...}







if __name__ == '__main__':


	##===========================
	##==== initialize dimensions
	##===========================
	#==== chromosome
	L = 100000000				# 100MB
	#==== individual
	n_individual = 185			# as in GTEx.v.4
	#==== tissue
	n_tissue = 17				# as in GTEx.v.4
	#==== SNP
	n_SNP = 100000				# 100K, so SNP density is 1/1,000bp
	#==== gene
	n_gene = 2000				# 2K, so gene density is 1/50,000bp
	#==== factor_cell
	n_factor_cell = 400			# as evaluated empirically
	#==== factor_batch (analogous to cell factor pathway)
	n_batch = 229
	n_batch_individual = 160		# as in GTEx.v.4
	n_batch_sample = 69			# as in GTEx.v.4
	n_factor_batch = 400			# as evaluated empirically, the same with cell factors





	##=======================================================
	##==== simulate variables (fill in the dimensions first)
	##=======================================================
	#==== chromosome
	#L = 100000000				# 100MB
	#==== individual
	#n_individual = 185			# as in GTEx.v.4
	#==== tissue
	#n_tissue = 17				# as in GTEx.v.4


	#==== SNP
	#n_SNP = 100000				# 100K, so SNP density is 1/1,000bp
	##SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
	for i in range(n_individual):
		SNP_rep[i] = []
		for j in range(n_SNP):
			SNP_rep[i].append(0)
	##SNP_pos_list = []
	for i in range(n_SNP):
		SNP_pos_list.append(0)
	##SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	# later on, after "SNP gene pos map"
	simu_genotype()


	#==== gene
	#n_gene = 2000				# 2K, so gene density is 1/50,000bp
	##gene_rep = {individual:{tissue:[], ...}, ...}			# gene expression list for tissue samples for individuals
	for i in range(n_individual):
		gene_rep[i] = {}
		for j in range(n_tissue):
			gene_rep[i][j] = []
			for k in range(n_gene):
				gene_rep[i][j].append(0)
	##gene_pos_list = []
	simu_gene()


	#==== SNP gene pos map
	#pos_map = {}							# {gene:[snp1, snp2], ...}
	SNP_gene_map()
	##SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	simu_genotype_beta()


	#==== factor_cell
	#n_factor_cell = 400			# as evaluated empirically
	##factor_cell_beta_rep = {cell factor:[], ...}			# coefficient lists for cell factors
	for i in range(n_factor_cell):
		factor_cell_beta_rep[i] = []
		for j in range(n_SNP):
			factor_cell_beta_rep[i].append(0)
		factor_cell_beta_rep[i].append(0)
	#beta_factor_cell_rep = {tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues
	for i in range(n_tissue):
		beta_factor_cell_rep[i] = {}
		for j in range(n_gene):
			beta_factor_cell_rep[i][j] = []
			for k in range(n_factor_cell):
				beta_factor_cell_rep[i][j].append(0)
			beta_factor_cell_rep[i][j].append(0)	# the constant item
	simu_cell_factor()


	#==== factor_batch (analogous to cell factor pathway)
	#n_batch = 229
	#n_batch_individual = 160		# as in GTEx.v.4
	#n_batch_sample = 69			# as in GTEx.v.4
	#n_factor_batch = 400
	##batch_individual_rep = {individual:[], ...}			# batch variable lists for individuals
	for i in range(n_individual):
		batch_individual_rep[i] = []
		for j in range(n_batch_individual):
			batch_individual_rep[i].append(0)
	##batch_tissue_sample_rep = {individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
	for i in range(n_individual):
		batch_tissue_sample_rep[i] = {}
		for j in range(n_tissue):
			batch_tissue_sample_rep[i][j] = []
			for k in range(n_batch_sample):
				batch_tissue_sample_rep[i][j].append(0)
	##factor_batch_beta_rep = {batch factor:[], ...}			# coefficient lists for batch factors
	for i in range(n_factor_batch):
		factor_batch_beta_rep[i] = []
		for j in range(n_batch):
			factor_batch_beta_rep[i].append(0)2
		factor_batch_beta_rep[i].append(0)
	##beta_factor_batch_rep = {gene:[], ...}
	for i in range(n_gene):
		beta_factor_batch_rep[i] = []
		for j in range(n_factor_batch):
			beta_factor_batch_rep[i].append(0)
		beta_factor_batch_rep[i].append(0)
	simu_batch_factor()












	##=======================================================
	##==== generate the expression profile through the model
	##=======================================================
	##gene_rep = {individual:{tissue:[], ...}, ...}			# gene expression list for tissue samples for individuals









	##======================
	##==== parameter saving
	##======================












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




