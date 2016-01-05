## Input: SNPs (independent sites), 0/1 binary array
## Output: expression level for some genes
## function: generate the coefficients according to the specified graphical model (the neural model in a Bayesian approach), and generate the corresponding expression level for ALL genes
## NOTE:
##	1. we should simulate tissue specificity and consistency (tissue hierarchy), as the modeling takes consideration of that;
##	2. when we say a factor, it's latent factor condensed from some original variables;
##	3. TODO: in next stage of simulation, we can use the true genotype with fake beta to generate the expression profile, as the genotype distribution (MAF) may not be exactly what we assume here;
##	4. we have all the coefficients one more item, which is the intercept of the regression: cis-SNP, SNP-cell factor, cell factor-gene, batch-batch factor, batch factor-gene
##	5. we might always simulate in the single chromosome mode
##	6. we might always simulate a full tensor, and then pick up the training set and testing set afterwards
##	7. the names of individuals/genes/tissues are indices coded in numbers from 0
##	8. there is no "gene_xymt.txt" information, so it's an empty file for the training program
##	x. ...




##================
##==== libraries
##================
import numpy as np
from utility import *
from simu_sub import *
from scipy.special import expit					# for logistic function



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



# DEBUG
#==== the gene expression signal from three pathways (cis, cell factor, batch effect)
gene_rep_cis = {}
gene_rep_cell = {}
gene_rep_batch = {}







if __name__ == '__main__':


	##===========================
	##==== initialize dimensions
	##===========================
	##==============================================================================================================
	# TODO: parameters tunable
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
	##==============================================================================================================






	##=============================================================
	##==== simulate variables (initialize space first if possible)
	##=============================================================

	# DEBUG
	print "now simulating genotype of individuals and their positions..."

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


	# DEBUG
	print "now simulating gene positions..."

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
	simu_gene_pos(gene_pos_list, L)


	# DEBUG
	print "now mapping the genes with SNPs, and simulating genotype beta..."

	#==== SNP gene pos map
	#pos_map = {}							# {gene:[snp1, snp2], ...}
	pos_map = {}
	SNP_gene_map(SNP_pos_list, gene_pos_list, pos_map)
	##SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	simu_genotype_beta(pos_map, n_tissue, n_gene, SNP_beta_rep)


	# DEBUG
	print "now simulating beta for SNP-cell factor pathway, and for cell factor-gene pathway..."

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


	# DEBUG
	print "now simulating batch variables and the beta of two pathways..."

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
		factor_batch_beta_rep[i].append(0)			# the constant item
	##beta_factor_batch_rep = {gene:[], ...}
	for i in range(n_gene):
		beta_factor_batch_rep[i] = []
		for j in range(n_factor_batch):
			beta_factor_batch_rep[i].append(0)
		beta_factor_batch_rep[i].append(0)			# the constant item
	simu_batch_factor(n_batch, n_batch_individual, n_batch_sample, n_factor_batch, batch_individual_rep, batch_tissue_sample_rep, factor_batch_beta_rep, beta_factor_batch_rep)




	# DEBUG
	print "now generating the expression level..."


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

	# DEBUG
	for i in range(n_individual):
		gene_rep_cis[i] = {}
		for j in range(n_tissue):
			gene_rep_cis[i][j] = []
			for k in range(n_gene):
				gene_rep_cis[i][j].append(0)
	for i in range(n_individual):
		gene_rep_cell[i] = {}
		for j in range(n_tissue):
			gene_rep_cell[i][j] = []
			for k in range(n_gene):
				gene_rep_cell[i][j].append(0)
	for i in range(n_individual):
		gene_rep_batch[i] = {}
		for j in range(n_tissue):
			gene_rep_batch[i][j] = []
			for k in range(n_gene):
				gene_rep_batch[i][j].append(0)


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
			#var_cell_hidden_factor[count] = value			# TODO: logistic twist; below
			var_cell_hidden_factor[count] = expit(value)

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
				#var_batch_hidden_factor[count] = value		# TODO: logistic twist; below
				var_batch_hidden_factor[count] = expit(value)

			for k in range(n_gene):

				# DEBUG
				# I will separately save the signal from the three pathway, to compare their intensity
				# the combined signal with noise added will be saved as usual

				y = 0
				#==== 1. cis- regulation
				y1 = 0
				index_start = pos_map[k][0]
				index_end = pos_map[k][1]
				for m in range(index_start, index_end + 1):
					dosage = SNP_rep[i][m]
					beta = SNP_beta_rep[j][k][m - index_start]
					y1 += dosage * beta
				y1 += 1 * SNP_beta_rep[j][k][-1]		# intercept
				gene_rep_cis[i][j][k] = y1
				y += y1

				#==== 2. cell factor regulation
				y2 = 0
				for count in range(n_factor_cell):
					factor = var_cell_hidden_factor[count]
					beta = beta_factor_cell_rep[j][k][count]
					y2 += factor * beta
				y2 += 1 * beta_factor_cell_rep[j][k][-1]	# intercept
				gene_rep_cell[i][j][k] = y2
				y += y2

				#==== 3. batch effect
				y3 = 0
				for count in range(n_factor_batch):
					factor = var_batch_hidden_factor[count]
					beta = beta_factor_batch_rep[k][count]
					y3 += factor * beta
				y3 += 1 * beta_factor_batch_rep[k][-1]		# intercept
				gene_rep_batch[i][j][k] = y3
				y += y3

				#==== 4. plus a small noise
				# Gaussian noise
				noise = np.random.normal(0, 1)
				y += noise

				#==== 5. assign this expression level to the gene
				gene_rep[i][j][k] = y









	# DEBUG
	print "now saving all the simulated data..."

	# NOTE: I will save two copies of the simulated data, one is in previous format, and another is in the format training program needs

	##=================================================
	##==== parameter saving (copy#1, the nature format)
	##=================================================
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

			#for rpkm in gene_rep[individual][tissue]:
			#	file_gene.write(str(rpkm) + "\t")
			# DEBUG
			for i in range(len(gene_rep[individual][tissue])):
				s = str(gene_rep[individual][tissue][i]) + " " + str(gene_rep_cis[individual][tissue][i]) + " " + str(gene_rep_cell[individual][tissue][i]) + " " + str(gene_rep_batch[individual][tissue][i])
				file_gene.write(s + "\t")

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
		file_batch_par_batch_hidden_gene.write(str(gene) + "\t")
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





	##=================================================================
	##==== parameter saving (copy#2, the format training program needs)
	##=================================================================
	# --> organized by tissue types (or individuals, for gene)
	file_list_individuals =			open("../simulation_data_reformat/list_individuals.txt", 'w')
	#file_SNP_var =				open("../simulation_data_reformat/genotype/chr1/SNP_dosage_IndividualID.txt", 'w')
	file_SNP_info =				open("../simulation_data_reformat/genotype/chr1/SNP_info.txt", 'w')
	#file_SNP_par = 			open("../simulation_data_reformat/SNP_par/TissueID.txt", 'w')
	file_batch_individuals = 		open("../simulation_data_reformat/batch_individuals.txt", 'w')
	file_batch_samples = 			open("../simulation_data_reformat/batch_samples.txt", 'w')
	file_batch_par_batch_batch_hidden = 	open("../simulation_data_reformat/batch_par_batch_batch_hidden.txt", 'w')
	file_batch_par_batch_hidden_gene = 	open("../simulation_data_reformat/batch_par_batch_hidden_gene.txt", 'w')
	file_cell_par_SNP_cell = 		open("../simulation_data_reformat/cell_par_SNP_cell.txt", 'w')
	#file_cell_par_cell_gene = 		open("../simulation_data_reformat/cell_par_cell_gene/TissueID.txt", 'w')
	file_expression = 			open("../simulation_data_reformat/expression.txt", 'w')
	file_gene_tss =				open("../simulation_data_reformat/gene_tss.txt", 'w')
	file_gene_SNP_map =			open("../simulation_data_reformat/gene_SNP_map.txt", 'w')


	# fill in the following: # -> need, ## -> doesn't
	#==== chromosome
	##L = 0								# unit: basepair

	#==== individual
	##n_individual = 0
	for i in range(n_individual):
		file_list_individuals.write(str(i) + '\n')

	#==== tissue
	##n_tissue = 0

	#==== SNP
	##n_SNP = 0
	#SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
	#SNP_pos_list = []
	#SNP_beta_rep = {} #{tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	for individual in SNP_rep:
		file_SNP_var =				open("../simulation_data_reformat/genotype/chr1/SNP_dosage_" + str(individual) + ".txt", 'w')
		for SNP in SNP_rep[individual]:
			file_SNP_var.write(str(SNP) + '\n')
		file_SNP_var.close()

	for i in range(len(SNP_pos_list)):
		SNP = str(i)
		pos = str(SNP_pos_list[i])
		file_SNP_info.write(SNP + ' ' + pos + '\n')

	for tissue in SNP_beta_rep:
		file_SNP_par = 				open("../simulation_data_reformat/SNP_par/" + str(tissue) + ".txt", 'w')
		for gene in SNP_beta_rep[tissue]:
			file_SNP_par.write(str(gene) + "\t")
			for beta in SNP_beta_rep[tissue][gene]:
				file_SNP_par.write(str(beta) + "\t")
			file_SNP_par.write("\n")
		file_SNP_par.close()

	#==== gene
	##n_gene = 0
	#gene_rep = {} #{individual:{tissue:[], ...}, ...}		# gene expression list for tissue samples for individuals
	#gene_pos_list = []

	# STEP (for expression.txt):
	# 1. get the individual-tissue list (this is the sample ID list) first
	# 2. get the expression matrix then
	# 3. save the two into the gene expression file
	# step#1:
	sample_list = []
	sample_index_map = {}
	count = 0
	for i in range(n_individual):
		for j in range(n_tissue):
			sample = str(i) + '-' + str(j)
			sample_list.append(sample)
			sample_index_map[sample] = count
			count += 1
	# step#2:
	n_sample = len(sample_list)
	expression_matrix = []
	for i in range(n_gene):
		expression_matrix.append([])
		for j in range(n_sample):
			expression_matrix[-1].append(0)
	for i in range(n_individual):
		for j in range(n_tissue):
			sample = str(i) + '-' + str(j)
			index = sample_index_map[sample]
			for k in range(n_gene):
				expression_matrix[k][index] = gene_rep[i][j][k]
	# step#3:
	file_expression.write("Name\tDescription\t")
	for sample in sample_list:
		file_expression.write(sample + '\t')
	file_expression.write('\n')
	for i in range(n_gene):
		file_expression.write(str(i) + '\t' + 'nan' + '\t')
		for j in range(len(expression_matrix[i])):
			rpkm = expression_matrix[i][j]
			file_expression.write(str(rpkm) + '\t')
		file_expression.write('\n')

	for i in range(len(gene_pos_list)):
		gene = str(i)
		chr = '1'		## NOTE: here we assume the single chromosome mode is always the case
		tss = str(gene_pos_list[i])
		file_gene_tss.write(gene + '\t' + chr + '\t' + tss + '\n')

	#==== SNP gene pos map
	#pos_map = {}							# {gene:[snp1, snp2], ...}
	for gene in pos_map:
		file_gene_SNP_map.write(str(gene) + "\t")
		file_gene_SNP_map.write(str(pos_map[gene][0]) + "\t")
		file_gene_SNP_map.write(str(pos_map[gene][1]))
		file_gene_SNP_map.write("\n")

	#==== factor_cell
	##n_factor_cell = 0
	##var_cell_hidden_factor = []					# n_factor_cell cell hidden factors to be filled in
	#factor_cell_beta_rep = {} #{cell factor:[], ...}		# coefficient lists for cell factors
	#beta_factor_cell_rep = {} #{tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues
	for factor in factor_cell_beta_rep:
		file_cell_par_SNP_cell.write(str(factor) + "\t")
		for beta in factor_cell_beta_rep[factor]:
			file_cell_par_SNP_cell.write(str(beta) + "\t")
		file_cell_par_SNP_cell.write("\n")
	for tissue in beta_factor_cell_rep:
		file_cell_par_cell_gene = 		open("../simulation_data_reformat/cell_par_cell_gene/" + str(tissue) + ".txt", 'w')
		for gene in beta_factor_cell_rep[tissue]:
			file_cell_par_cell_gene.write(str(gene) + "\t")
			for beta in beta_factor_cell_rep[tissue][gene]:
				file_cell_par_cell_gene.write(str(beta) + "\t")
			file_cell_par_cell_gene.write("\n")
		file_cell_par_cell_gene.close()

	#==== factor_batch (analogous to cell factor pathway)
	##n_batch = 0							# n_batch = n_batch_individual + n_batch_sample
	##n_batch_individual = 0
	##n_batch_sample = 0
	##n_factor_batch = 0						# the number of latent batch factors
	##var_batch_hidden_factor = []					# n_factor_batch batch hidden factors to be filled in
	#batch_individual_rep = {} #{individual:[], ...}			# batch variable lists for individuals
	#batch_tissue_sample_rep = {} #{individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
	#factor_batch_beta_rep = {} #{batch factor:[], ...}		# coefficient lists for batch factors
	#beta_factor_batch_rep = {} #{gene:[], ...}			# coefficient lists of batch factors for genes
	file_batch_individuals.write("individual\t")
	for i in range(n_batch_individual):
		file_batch_individuals.write(str(i) + "\t")
	file_batch_individuals.write("\n")
	for individual in batch_individual_rep:
		file_batch_individuals.write(str(individual) + "\t")
		for value in batch_individual_rep[individual]:
			file_batch_individuals.write(str(value) + "\t")
		file_batch_individuals.write("\n")

	file_batch_samples.write("sample\t")
	for i in range(n_batch_sample):
		file_batch_samples.write(str(i) + "\t")
	file_batch_samples.write("\n")
	for individual in batch_tissue_sample_rep:			# a sample is defined for per individual per tissue
		for tissue in batch_tissue_sample_rep[individual]:
			sample = str(individual) + "-" + str(tissue)
			file_batch_samples.write(sample + "\t")
			for value in batch_tissue_sample_rep[individual][tissue]:
				file_batch_samples.write(str(value) + "\t")
			file_batch_samples.write("\n")

	for factor in factor_batch_beta_rep:
		file_batch_par_batch_batch_hidden.write(str(factor) + "\t")
		for beta in factor_batch_beta_rep[factor]:
			file_batch_par_batch_batch_hidden.write(str(beta) + "\t")
		file_batch_par_batch_batch_hidden.write("\n")

	for gene in beta_factor_batch_rep:
		file_batch_par_batch_hidden_gene.write(str(gene) + "\t")
		for beta in beta_factor_batch_rep[gene]:
			file_batch_par_batch_hidden_gene.write(str(beta) + "\t")
		file_batch_par_batch_hidden_gene.write("\n")





	file_list_individuals.close()# =			open("../simulation_data_reformat/list_individuals.txt", 'w')
	#file_SNP_var =						open("../simulation_data_reformat/genotype/chr1/SNP_dosage_IndividualID.txt", 'w')
	file_SNP_info.close()# =				open("../simulation_data_reformat/genotype/chr1/SNP_info.txt", 'w')
	#file_SNP_par = 					open("../simulation_data_reformat/SNP_par/TissueID.txt", 'w')
	file_batch_individuals.close()# = 			open("../simulation_data_reformat/batch_individuals.txt", 'w')
	file_batch_samples.close()# = 				open("../simulation_data_reformat/batch_samples.txt", 'w')
	file_batch_par_batch_batch_hidden.close()# = 		open("../simulation_data_reformat/batch_par_batch_batch_hidden.txt", 'w')
	file_batch_par_batch_hidden_gene.close()# = 		open("../simulation_data_reformat/batch_par_batch_hidden_gene.txt", 'w')
	file_cell_par_SNP_cell.close()# = 			open("../simulation_data_reformat/cell_par_SNP_cell.txt", 'w')
	#file_cell_par_cell_gene = 				open("../simulation_data_reformat/cell_par_cell_gene/TissueID.txt", 'w')
	file_expression.close()# = 				open("../simulation_data_reformat/expression.txt", 'w')
	file_gene_tss.close()# =				open("../simulation_data_reformat/gene_tss.txt", 'w')
	file_gene_SNP_map.close()# =				open("../simulation_data_reformat/gene_SNP_map.txt", 'w')



