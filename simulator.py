## Input: SNPs (independent sites), 0/1 binary array
## Output: expression level for some genes
## function: generate the coefficients according to the specified graphical model (the neural model in a Bayesian approach), and generate the corresponding expression level for ALL genes
## notes:
##	1. we should simulate tissue specificity, as the modeling takes consideration of that;
##	2. when we say a factor, it's latent factor condensed from some original variables;
##	3. TODO: in next stage of simulation, we can use the true genotype with fake beta to generate the expression profile, as the genotype distribution (MAF) may not be exactly what we assume here;
##	4. xxx



##================
##==== libraries
##================
import numpy as np
from utility import *
from simu_sub import *







##=====================
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
SNP_beta_rep = {} #{tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess

#==== gene
n_gene = 0
gene_rep = {} #{individual:{tissue:[], ...}, ...}		# gene expression list for tissue samples for individuals
gene_pos_list = []

#==== SNP gene pos map
pos_map = {}							# {gene:[snp1, snp2], ...}

#==== factor_cell
n_factor_cell = 0
factor_cell_beta_rep = {} #{cell factor:[], ...}		# coefficient lists for cell factors
beta_factor_cell_rep = {} #{tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues

#==== factor_batch (analogous to cell factor pathway)
n_batch = 0							# n_batch = n_batch_individual + n_batch_sample
n_batch_individual = 0
n_batch_sample = 0
n_factor_batch = 0						# the number of latent batch factors
batch_individual_rep = {} #{individual:[], ...}			# batch variable lists for individuals
batch_tissue_sample_rep = {} #{individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
factor_batch_beta_rep = {} #{batch factor:[], ...}		# coefficient lists for batch factors
beta_factor_batch_rep = {} #{gene:[], ...}			# coefficient lists of batch factors for genes








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
	n_factor_batch = 50			# as evaluated empirically





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
	simu_gene_pos(gene_pos_list)


	#==== SNP gene pos map
	#pos_map = {}							# {gene:[snp1, snp2], ...}
	SNP_gene_map(SNP_pos_list, gene_pos_list, pos_map)
	##SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	simu_genotype_beta(pos_map, n_tissue, n_gene, SNP_beta_rep)



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
	simu_cell_factor(factor_cell_beta_rep, n_SNP, beta_factor_cell_rep, n_factor_cell)



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
	simu_batch_factor(n_batch, n_batch_individual, n_batch_sample, n_factor_batch, batch_individual_rep, batch_tissue_sample_rep, factor_batch_beta_rep, beta_factor_batch_rep)








	##=======================================================
	##==== generate the expression profile through the model
	##=======================================================
	##gene_rep = {individual:{tissue:[], ...}, ...}			# gene expression list for tissue samples for individuals
	for i in range(n_individual):
		gene_rep[i] = {}
		for j in range(n_tissue):
			gene_rep[i][j] = []
			for k in range(n_gene):
				gene_rep[i][j].append(0)

	# three pathway: 1. cis- regulation; 2. trans- regulation; 3. batch effect
	for i in range(n_individual):		#
		for j in range(n_tissue):	#
			for k in range(n_gene):	#
				y = 0
				#==== 1. cis- regulation
				index_start = pos_map[k][0]
				index_end = pos_map[k][1]
				for m in range(index_start, index_end+1):
					dosage = SNP_rep[i][m]
					beta = SNP_beta_rep[j][k][m]
					y += dosage * beta
				y += SNP_beta_rep[j][k][-1]	# intercept

				#==== 2. trans- regulation



				#==== 3. batch effect

				

				noise = np.random.normal()
				y += noise
				gene_rep[i][j][k] = y




	"""
	#==== individual
	n_individual = 0
	#==== tissue
	n_tissue = 0
	#==== SNP
	n_SNP = 0
	SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
	SNP_beta_rep = {} #{tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	#==== gene
	n_gene = 0
	gene_rep = {} #{individual:{tissue:[], ...}, ...}		# gene expression list for tissue samples for individuals
	gene_pos_list = []
	#==== SNP gene pos map
	pos_map = {}							# {gene:[snp1, snp2], ...}
	#==== factor_cell
	n_factor_cell = 0
	factor_cell_beta_rep = {} #{cell factor:[], ...}		# coefficient lists for cell factors
	beta_factor_cell_rep = {} #{tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues
	#==== factor_batch (analogous to cell factor pathway)
	n_batch = 0							# n_batch = n_batch_individual + n_batch_sample
	n_batch_individual = 0
	n_batch_sample = 0
	n_factor_batch = 0						# the number of latent batch factors
	batch_individual_rep = {} #{individual:[], ...}			# batch variable lists for individuals
	batch_tissue_sample_rep = {} #{individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
	factor_batch_beta_rep = {} #{batch factor:[], ...}		# coefficient lists for batch factors
	beta_factor_batch_rep = {} #{gene:[], ...}			# coefficient lists of batch factors for genes
	"""














	##======================
	##==== parameter saving
	##======================














