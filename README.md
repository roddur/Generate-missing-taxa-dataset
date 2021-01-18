# Generate-missing-taxa-dataset

## Description
Deletes clades from gene trees based on $M_{clade}$ model described in -
Nute, Michael, et al. "The performance of coalescent-based species tree estimation methods under models of missing data." BMC genomics 19.5 (2018): 1-22.

## Detailed Description
First, it lists all the clades from the species tree having the number of species above and below two thresholds. Next, the gene trees from which a clade cannot be deleted (have enough leaves and a clade has been found contained entirely in the gene tree after a specified number of tries) are filtered out.From the remaining gene trees, a certain percentage of all trees have been chosen and from each tree, a random clade listed from the species tree has been deleted. The process is the same as described as Nute et. al. when all the gene trees have the same number of species.
