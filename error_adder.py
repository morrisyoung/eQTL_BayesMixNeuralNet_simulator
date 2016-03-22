## func: add some errors to the generated models (parameters)
#
## Feb.4: it seems not necessary to do this, as we are going wild even we start from the real parameter values
#
## Feb.8: we will do this, to see how the likelihood curve changes; we should do this anyway
#
## notes:
#	1. we should run this after the simulator routine;



##================
import numpy as np
import os




if __name__ == "__main__":


	##======== create the folder first of all
	if not os.path.isdir("../simulation_para_with_error"):
		os.mkdir("../simulation_para_with_error")
		os.mkdir("../simulation_para_with_error/para_init_cis_gene")
		os.mkdir("../simulation_para_with_error/para_init_cellenv_gene")
	else:
		if not os.path.isdir("../simulation_para_with_error/para_init_cis_gene"):
			os.mkdir("../simulation_para_with_error/para_init_cis_gene")
		if not os.path.isdir("../simulation_para_with_error/para_init_cellenv_gene"):
			os.mkdir("../simulation_para_with_error/para_init_cellenv_gene")


	##========= tissue file mapping
	tissue_file_mapping = {}
	file = open("../simulation_para/etissue_list_init.txt", 'r')
	file1 = open("../simulation_para_with_error/etissue_list_init.txt", 'w')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split("\t")
		tissue = line[0]
		index = line[1]

		tissue_file_mapping[tissue] = index
		file1.write(tissue + '\t' + index + '\n')

	file.close()
	file1.close()




	##========= cis_gene
	# TODO: parameter tunable
	lamb = 0
	std = 1

	for tissue in tissue_file_mapping:
		tissue_index = tissue_file_mapping[tissue]

		filename = "etissue" + tissue_index + ".txt"
		file = open("../simulation_para/para_init_cis_gene/" + filename, 'r')
		file1 = open("../simulation_para_with_error/para_init_cis_gene/" + filename, 'w')

		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			gene = line[0]			# NOTE: why we need the gene ID to index its parameters?
			file1.write(gene + '\t')

			line = map(lambda x: float(x), line[1:])
			for beta in line:
				# add error
				beta += np.random.normal(lamb, std)
				file1.write(str(beta) + '\t')
			file1.write('\n')

		file.close()
		file1.close()



	##========= cellenv_gene
	# TODO: parameter tunable
	lamb = 0
	std = 1

	for tissue in tissue_file_mapping:
		tissue_index = tissue_file_mapping[tissue]

		filename = "etissue" + tissue_index + ".txt"
		file = open("../simulation_para/para_init_cellenv_gene/" + filename, 'r')
		file1 = open("../simulation_para_with_error/para_init_cellenv_gene/" + filename, 'w')

		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			line = map(lambda x: float(x), line)
			for beta in line:
				# add error
				beta += np.random.normal(lamb, std)
				file1.write(str(beta) + '\t')
			file1.write('\n')

		file.close()
		file1.close()



	##========= snp_cellenv
	# TODO: parameter tunable
	lamb = 0
	std = 1

	file = open("../simulation_para/para_init_snp_cellenv.txt", 'r')
	file1 = open("../simulation_para_with_error/para_init_snp_cellenv.txt", 'w')

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		line = map(lambda x: float(x), line)
		for beta in line:
			# add error
			beta += np.random.normal(lamb, std)
			file1.write(str(beta) + '\t')
		file1.write('\n')

	file.close()
	file1.close()



	##========= batch_batch_hidden
	# TODO: parameter tunable
	lamb = 0
	std = 1

	file = open("../simulation_para/para_init_batch_batch_hidden.txt", 'r')
	file1 = open("../simulation_para_with_error/para_init_batch_batch_hidden.txt", 'w')

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		line = map(lambda x: float(x), line)
		for beta in line:
			# add error
			beta += np.random.normal(lamb, std)
			file1.write(str(beta) + '\t')
		file1.write('\n')

	file.close()
	file1.close()



	##========= batch_hidden_gene
	# TODO: parameter tunable
	lamb = 0
	std = 1

	file = open("../simulation_para/para_init_batch_hidden_gene.txt", 'r')
	file1 = open("../simulation_para_with_error/para_init_batch_hidden_gene.txt", 'w')

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		line = map(lambda x: float(x), line)
		for beta in line:
			# add error
			beta += np.random.normal(lamb, std)
			file1.write(str(beta) + '\t')
		file1.write('\n')

	file.close()
	file1.close()




	print "error adding done."


