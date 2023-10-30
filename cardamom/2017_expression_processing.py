import numpy as np
import pandas as pd
import copy

# Get array of all genes of interest by their names
interest_genes = np.loadtxt("dicty/Data/2017/jump_genes.csv", delimiter = ',', dtype = str, encoding="utf_8_sig")

# Remove genes beginning with DDB
interest_genes = interest_genes[~np.char.startswith(interest_genes, 'DDB')]

# Get array of gene names corresponding to gene IDs
ids = pd.read_csv("dicty/Data/updated_gene_IDs.csv", header=None)

# Create dictionary mapping gene names to gene IDs
gene_ids = []
gene_name_to_id = pd.Series(ids[0].values, index=ids[1]).to_dict()

# Get array of all gene IDs of interest
for i in range(len(interest_genes)):
    gene_ids.append(gene_name_to_id[interest_genes[i]])
gene_ids = np.array(gene_ids)

# Get expression matrix for all genes
full_data = np.loadtxt("dicty/Data/2017/norm_data.csv", delimiter = ',', dtype = str)

# Get row indices of gene of interests in expression matrix
genes = copy.deepcopy(full_data[:, 0])
gene_indices = [0] + [np.where(genes==item)[0][0] for item in gene_ids]



# Get expression matrix for genes of interest
full_data[0, 0] = 0
data = full_data[gene_indices]

# Replace first column with indices from 0 to N
data[1:, 0] = np.arange(1, len(data))
data=data.astype(float)
data=data.astype(int)

# Insert stimulus into expression matrix and gene list
data = np.insert(data,1, np.zeros((len(data[0]))), axis=0)
for i in range(len(data[0])):
    if data[0, i] != 0:
        data[1, i] = 1

# Add gene id index
interest_genes = np.insert(interest_genes.astype(object), 0, "Stimulus")
interest_genes = np.column_stack((np.arange(len(interest_genes)), interest_genes))

np.savetxt("dicty/Data/panel_real.txt", data, delimiter = '\t', fmt = '%i')
np.savetxt("dicty/Data/panel_genes.txt", interest_genes, delimiter = '\t', fmt = '%s')

