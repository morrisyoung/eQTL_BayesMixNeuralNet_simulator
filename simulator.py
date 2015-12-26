## Input: SNPs (independent sites), 0/1 binary array
## Output: expression level for some genes
## function: generate the coefficients according to the specified graphical model (the neural model in a Bayesian approach), and generate the corresponding expression level for ALL genes
## notes:
##	1. we should simulate tissue specificity, as the modeling takes consideration of that;
##	2. when we say a factor, it's latent factor condensed from some original variables;
##	3. TODO: in next stage of simulation, we can use the true genotype with fake beta to generate the expression profile, as the genotype distribution (MAF) may not be exactly what we assume here;
##	4. we should simulate the tissue hierarchy, as we need to integrate that into the framework (learning)
##	5. xxx




##================
##==== libraries
##================
import numpy as np
from utility import *
from simu_sub import *



##=====================
##==== global variables
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
var_cell_hidden_factor = []					# n_factor_cell cell hidden factors to be filled in
factor_cell_beta_rep = {} #{cell factor:[], ...}		# coefficient lists for cell factors
beta_factor_cell_rep = {} #{tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues

#==== factor_batch (analogous to cell factor pathway)
n_batch = 0							# n_batch = n_batch_individual + n_batch_sample
n_batch_individual = 0
n_batch_sample = 0
n_factor_batch = 0						# the number of latent batch factors
var_batch_hidden_factor = []					# n_factor_batch batch hidden factors to be filled in
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




	##=============================================================
	##==== simulate variables (initialize space first if possible)
	##=============================================================
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
	# later on, after "SNP gene pos map", as we now don't know the mapping of cis- SNPs
	simu_genotype(SNP_rep, SNP_pos_list, L)


	#==== gene
	#n_gene = 2000				# 2K, so gene density is 1/50,000bp
	##gene_rep = {individual:{tissue:[], ...}, ...}			# gene expression list for tissue samples for individuals
	#for i in range(n_individual):
	#	gene_rep[i] = {}
	#	for j in range(n_tissue):
	#		gene_rep[i][j] = []
	#		for k in range(n_gene):
	#			gene_rep[i][j].append(0)
	##gene_pos_list = []
	for i in range(n_gene):
		gene_pos_list.append(0)
	simu_gene_pos(gene_pos_list)


	#==== SNP gene pos map
	#pos_map = {}							# {gene:[snp1, snp2], ...}
	SNP_gene_map(SNP_pos_list, gene_pos_list, pos_map)
	##SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	simu_genotype_beta(pos_map, n_tissue, n_gene, SNP_beta_rep)


	#==== factor_cell
	#n_factor_cell = 400			# as evaluated empirically
	##var_cell_hidden_factor = []					# n_factor_cell cell hidden factors to be filled in
	##factor_cell_beta_rep = {cell factor:[], ...}			# coefficient lists for cell factors
	for i in range(n_factor_cell):
		factor_cell_beta_rep[i] = []
		for j in range(n_SNP):
			factor_cell_beta_rep[i].append(0)
		factor_cell_beta_rep[i].append(0)			# the constant item
	#beta_factor_cell_rep = {tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues
	for i in range(n_tissue):
		beta_factor_cell_rep[i] = {}
		for j in range(n_gene):
			beta_factor_cell_rep[i][j] = []
			for k in range(n_factor_cell):
				beta_factor_cell_rep[i][j].append(0)
			beta_factor_cell_rep[i][j].append(0)		# the constant item
	simu_cell_factor(factor_cell_beta_rep, n_SNP, beta_factor_cell_rep, n_factor_cell)


	#==== factor_batch (analogous to cell factor pathway)
	#n_batch = 229
	#n_batch_individual = 160		# as in GTEx.v.4
	#n_batch_sample = 69			# as in GTEx.v.4
	#n_factor_batch = 400
	##var_batch_hidden_factor = []					# n_factor_batch batch hidden factors to be filled in
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
			factor_batch_beta_rep[i].append(0)
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

	##var_cell_hidden_factor = []					# n_factor_cell cell hidden factors to be filled in
	for i in range(n_factor_cell):
		var_cell_hidden_factor.append(0)

	##var_batch_hidden_factor = []					# n_factor_batch batch hidden factors to be filled in
	for i in range(n_factor_batch):
		var_batch_hidden_factor.append(0)

	# three pathway: 1. cis- regulation; 2. trans- regulation; 3. batch effect
	for i in range(n_individual):
		# get the cell hidden factors for this individual
		for count in range(n_factor_cell):
			value = 0
			for count1 in range(n_SNP):
				value += SNP_rep[i][count1] * factor_cell_beta_rep[count][count1]
			value += 1 * factor_cell_beta_rep[count][-1]
			var_cell_hidden_factor[count] = value

		for j in range(n_tissue):
			# get the batch hidden factors for this individual in this tissue
			batch_list = []
			for count in range(n_batch_individual):
				batch_list.append(batch_individual_rep[i][count])
			for count in range(n_batch_sample):
				batch_list.append(batch_tissue_sample_rep[i][j][count])
			for count in range(n_factor_batch):
				value = 0
				for count1 in range(n_batch):
					value += batch_list[count1] * factor_batch_beta_rep[count][count1]
				value += 1 * factor_batch_beta_rep[count][-1]
				var_batch_hidden_factor[count] = value

			for k in range(n_gene):
				y = 0
				#==== 1. cis- regulation
				index_start = pos_map[k][0]
				index_end = pos_map[k][1]
				for m in range(index_start, index_end + 1):
					dosage = SNP_rep[i][m]
					beta = SNP_beta_rep[j][k][m - index_start]
					y += dosage * beta
				y += 1 * SNP_beta_rep[j][k][-1]			# intercept

				#==== 2. cell factor regulation
				for count in range(n_factor_cell):
					factor = var_cell_hidden_factor[count]
					beta = beta_factor_cell_rep[j][k][count]
					y += factor * beta
				y += 1 * beta_factor_cell_rep[j][k][-1]		# intercept

				#==== 3. batch effect
				for count in range(n_factor_batch):
					factor = var_batch_hidden_factor[count]
					beta = beta_factor_batch_rep[k][count]
					y += factor * beta
				y += 1 * beta_factor_batch_rep[k][-1]		# intercept

				#==== 4. plus a small noise
				# Gaussian noise
				noise = np.random.normal(0, 1)
				y += noise

				#==== 5. assign this expression level to the gene
				gene_rep[i][j][k] = y





	##======================
	##==== parameter saving
	##======================
	# --> organized by tissue types (or individuals, for gene)
	file_meta =				open("../simulation_data/meta.txt", 'w')
	file_SNP_var =				open("../simulation_data/SNP_var.txt", 'w')
	#file_SNP_par = 				open("../simulation_data/SNP_par/xxx.txt", 'w')
	file_batch_var_individual = 		open("../simulation_data/batch_var_individual.txt", 'w')
	file_batch_var_sample = 		open("../simulation_data/batch_var_sample.txt", 'w')
	file_batch_par_batch_batch_hidden = 	open("../simulation_data/batch_par_batch_batch_hidden.txt", 'w')
	file_batch_par_batch_hidden_gene = 	open("../simulation_data/batch_par_batch_hidden_gene.txt", 'w')
	file_cell_par_SNP_cell = 		open("../simulation_data/cell_par_SNP_cell.txt", 'w')
	#file_cell_par_cell_gene = 		open("../simulation_data/cell_par_cell_gene/xxx.txt", 'w')
	#file_gene = 				open("../simulation_data/gene/xxx.txt", 'w')
	file_gene_SNP_map =			open("../simulation_data/gene_SNP_map.txt", 'w')


	# fill in the following: # -> need, ## -> doesn't
	#==== chromosome
	#L = 0								# unit: basepair
	file_meta.write("chromosome\t" + str(L) + "\n")

	#==== individual
	#n_individual = 0
	file_meta.write("n_individual\t" + str(n_individual) + "\n")

	#==== tissue
	#n_tissue = 0
	file_meta.write("n_tissue\t" + str(n_tissue) + "\n")

	#==== SNP
	#n_SNP = 0
	#SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
	#SNP_pos_list = []
	#SNP_beta_rep = {} #{tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	file_meta.write("n_SNP\t" + str(n_SNP) + "\n")
	for individual in SNP_rep:
		file_SNP_var.write(str(individual) + "\t")
		for SNP in SNP_rep[individual]:
			file_SNP_var.write(str(SNP) + "\t")
		file_SNP_var.write("\n")
	file_meta.write("SNP_pos_list\t")
	for pos in SNP_pos_list:
		file_meta.write(str(pos) + "\t")
	file_meta.write("\n")
	for tissue in SNP_beta_rep:
		file_SNP_par = 				open("../simulation_data/SNP_par/" + str(tissue) + ".txt", 'w')
		for gene in SNP_beta_rep[tissue]:
			file_SNP_par.write(str(gene) + "\t")
			for beta in SNP_beta_rep[tissue][gene]:
				file_SNP_par.write(str(beta) + "\t")
			file_SNP_par.write("\n")
		file_SNP_par.close()

	#==== gene
	#n_gene = 0
	#gene_rep = {} #{individual:{tissue:[], ...}, ...}		# gene expression list for tissue samples for individuals
	#gene_pos_list = []
	file_meta.write("n_gene\t" + str(n_gene) + "\n")
	for individual in gene_rep:
		file_gene = 				open("../simulation_data/gene/" + str(individual) + ".txt", 'w')
		for tissue in gene_rep[individual]:
			file_gene.write(str(tissue) + "\t")
			for rpkm in gene_rep[individual][tissue]:
				file_gene.write(str(rpkm) + "\t")
			file_gene.write("\n")
		file_gene.close()
	file_meta.write("gene_pos_list\t")
	for pos in gene_pos_list:
		file_meta.write(str(pos) + "\t")
	file_meta.write("\n")

	#==== SNP gene pos map
	#pos_map = {}							# {gene:[snp1, snp2], ...}
	for gene in pos_map:
		file_gene_SNP_map.write(str(gene) + "\t")
		file_gene_SNP_map.write(str(pos_map[gene][0]) + "\t")
		file_gene_SNP_map.write(str(pos_map[gene][1]))
		file_gene_SNP_map.write("\n")

	#==== factor_cell
	#n_factor_cell = 0
	##var_cell_hidden_factor = []					# n_factor_cell cell hidden factors to be filled in
	#factor_cell_beta_rep = {} #{cell factor:[], ...}		# coefficient lists for cell factors
	#beta_factor_cell_rep = {} #{tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues
	file_meta.write("n_factor_cell\t" + str(n_factor_cell) + "\n")
	for factor in factor_cell_beta_rep:
		file_cell_par_SNP_cell.write(str(factor) + "\t")
		for beta in factor_cell_beta_rep[factor]:
			file_cell_par_SNP_cell.write(str(beta) + "\t")
		file_cell_par_SNP_cell.write("\n")
	for tissue in beta_factor_cell_rep:
		file_cell_par_cell_gene = 		open("../simulation_data/cell_par_cell_gene/" + str(tissue) + ".txt", 'w')
		for gene in beta_factor_cell_rep[tissue]:
			file_cell_par_cell_gene.write(str(gene) + "\t")
			for beta in beta_factor_cell_rep[tissue][gene]:
				file_cell_par_cell_gene.write(str(beta) + "\t")
			file_cell_par_cell_gene.write("\n")
		file_cell_par_cell_gene.close()

	#==== factor_batch (analogous to cell factor pathway)
	#n_batch = 0							# n_batch = n_batch_individual + n_batch_sample
	#n_batch_individual = 0
	#n_batch_sample = 0
	#n_factor_batch = 0						# the number of latent batch factors
	##var_batch_hidden_factor = []					# n_factor_batch batch hidden factors to be filled in
	#batch_individual_rep = {} #{individual:[], ...}			# batch variable lists for individuals
	#batch_tissue_sample_rep = {} #{individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
	#factor_batch_beta_rep = {} #{batch factor:[], ...}		# coefficient lists for batch factors
	#beta_factor_batch_rep = {} #{gene:[], ...}			# coefficient lists of batch factors for genes
	file_meta.write("n_batch\t" + str(n_batch) + "\n")
	file_meta.write("n_batch_individual\t" + str(n_batch_individual) + "\n")
	file_meta.write("n_batch_sample\t" + str(n_batch_sample) + "\n")
	file_meta.write("n_factor_batch\t" + str(n_factor_batch) + "\n")
	for individual in batch_individual_rep:
		file_batch_var_individual.write(str(individual) + "\t")
		for value in batch_individual_rep[individual]:
			file_batch_var_individual.write(str(value) + "\t")
		file_batch_var_individual.write("\n")
	for individual in batch_tissue_sample_rep:			# a sample is defined for per individual per tissue
		for tissue in batch_tissue_sample_rep[individual]:
			key = str(individual) + "-" + str(tissue)
			file_batch_var_sample.write(key + "\t")
			for value in batch_tissue_sample_rep[individual][tissue]:
				file_batch_var_sample.write(str(value) + "\t")
			file_batch_var_sample.write("\n")
	for factor in factor_batch_beta_rep:
		file_batch_par_batch_batch_hidden.write(str(factor) + "\t")
		for beta in factor_batch_beta_rep[factor]:
			file_batch_par_batch_batch_hidden.write(str(beta) + "\t")
		file_batch_par_batch_batch_hidden.write("\n")
	for gene in beta_factor_batch_rep:
		file_batch_par_batch_hidden_gene.write(str(factor) + "\t")
		for beta in beta_factor_batch_rep[gene]:
			file_batch_par_batch_hidden_gene.write(str(beta) + "\t")
		file_batch_par_batch_hidden_gene.write("\n")



	file_meta.close()# =				open("../simulation_data/meta.txt", 'w')
	file_SNP_var.close()# =				open("../simulation_data/SNP_var.txt", 'w')
	#file_SNP_par = 				open("../simulation_data/SNP_par/xxx.txt", 'w')
	file_batch_var_individual.close()# = 		open("../simulation_data/batch_var_individual.txt", 'w')
	file_batch_var_sample.close()# = 		open("../simulation_data/batch_var_sample.txt", 'w')
	file_batch_par_batch_batch_hidden.close()# = 	open("../simulation_data/batch_par_batch_batch_hidden.txt", 'w')
	file_batch_par_batch_hidden_gene.close()# = 	open("../simulation_data/batch_par_batch_hidden_gene.txt", 'w')
	file_cell_par_SNP_cell.close()# = 		open("../simulation_data/cell_par_SNP_cell.txt", 'w')
	#file_cell_par_cell_gene = 		open("../simulation_data/cell_par_cell_gene/xxx.txt", 'w')
	#file_gene = 				open("../simulation_data/gene/xxx.txt", 'w')
	file_gene_SNP_map.close()# =			open("../simulation_data/gene_SNP_map.txt", 'w')

