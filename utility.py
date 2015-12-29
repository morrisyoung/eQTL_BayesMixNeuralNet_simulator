## functions in this source code:
##	1. map the cis- regions of all genes;
##	2. xxx





## input: the pos of all SNPs; the pos of all genes
## output: the cis- region index (start, end) of SNPs for all genes
def SNP_gene_map(SNP_pos_list, gene_pos_list, pos_map):

	for i in range(len(gene_pos_list)):
		index_gene = i
		index_start = 0
		index_end = 0

		pos_gene = gene_pos_list[index_gene]
		index = 0
		while (SNP_pos_list[index] - pos_gene) < -1000000:
			index += 1
		index_start = index
		index += 1
		while (SNP_pos_list[index] - pos_gene < 1000000):
			index += 1
			if index == len(SNP_pos_list):
				break
		index -= 1
		index_end = index
		pos_map[index_gene] = [index_start, index_end]

	return


