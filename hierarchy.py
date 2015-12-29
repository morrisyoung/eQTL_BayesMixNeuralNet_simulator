## functions: all the components for tissue hierarchy regulation coefficients generation

import numpy as np


hierarchy = {}



def HierarchyInit():
	global hierarchy
	# initialize the hierarchy (graph structure, and branch length)
	# let's design how the hierarchy is stored first:
	#	1. we hash root and all other internal nodes and leaves, all the internal nodes (including the root) point to their two childen, while all the leaves point to NULL;
	#	2. we traverse the hierarchy and simulate the childen accordingly;
	#	3. indices of nodes are the hashing keys, and hashing values include: [value array, I[is root], I[is leaf], indices of two childen, branch lengths to two childen, index of parent], which is also the file format to save this hierarchy;
	#	4. TODO: we need another index file to map indices of leaves and tissue types, for training program use

	#==== load the hierarchy
	file = open("./hierarchy_tissue.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		key = int(line[0])
		#values = line[1]
		root = int(line[2])
		leaf = int(line[3])
		childen = map(lambda x: int(x), line[4].split(' '))
		branch = map(lambda x: float(x), line[5].split(' '))
		parent = int(line[6])

		hierarchy[key] = [[], root, leaf, childen, branch, parent]
	file.close()

	return



def Random_Gaussian(mean, var):
	# generate random number from specified Gaussian
	mu = mean
	sigma = np.sqrt(var)
	number = np.random.normal(mu, sigma, 1)
	return number



## hierarchy:
##	key: [value array, I[is root], I[is leaf], indices of two childen, branch lengths to two childen, index of parent]
def Hierarchy_NodeFill(node, length):
	global hierarchy

	#==== generate the current value array
	# root is from 0 mean unit variance
	if hierarchy[node][1]:
		list_temp = [0] * length
		for count in range(len(list_temp)):
			number = Random_Gaussian(0, 1)
			list_temp[count] = number

		hierarchy[node][0] = list_temp
	else:
		parent = hierarchy[node][5]
		mean = hierarchy[parent][0]			# array
		var = 0						# scaler
		if node == hierarchy[parent][3][0]:		# the left child of its parent
			var = hierarchy[parent][4][0]
		else:						# the right child of its parent
			var = hierarchy[parent][4][1]

		list_temp = [0] * length
		for count in range(len(list_temp)):
			number = Random_Gaussian(mean[count], var)
			list_temp[count] = number

		hierarchy[node][0] = list_temp

	#==== traverse the childen if there are
	if not hierarchy[node][2]:
		# left child
		child = hierarchy[node][3][0]
		Hierarchy_NodeFill(child, length)
		# right child
		child = hierarchy[node][3][1]
		Hierarchy_NodeFill(child, length)

	#==== return either to the upper stack, or the main program
	return
	



## target:
##	output all the leaves (with given array length) of a given hierarchy (structure, and branch length), starting from 0 mean and unit variance
## notes:
##	1. here we will use the uni-variance Gaussian (assuming all the regulators are independent; we can guarantee this for cis- SNPs, as we have pruned some correlated SNPs);
##	2. so the final MAP includes two parts: a. the gradient descent (or Gibbs sampling) for the regression part, b. the weighted least square for the hierarchy part
##	3. in this simulation, we should just simulate the leaves according to the weighted hierarchy
##	4. xxx
def HierarchyCal(length):
	global hierarchy
	leaves = {}	# there are #tissue elements inside, representing the regulation coefficients from each tissue
	
	#==== update the hierarchy (mainly for the value array item)
	root = 0
	for node in hierarchy:
		if not hierarchy[node][1]:
			continue
		else:
			root = node
			break

	Hierarchy_NodeFill(root, length)

	#==== get out all the leaves, drop in "leaves", and return it
	for node in hierarchy:
		if hierarchy[node][2]:
			leaves[node] = hierarchy[node][0]

	return leaves


