import numpy as np
import pandas as pd

# Get all degradation rates
rates = np.loadtxt("dicty/Data/degradation_rates.csv", delimiter = ',', dtype = str, encoding="utf_8_sig")
rates = rates[1:, [0, 12]]


# Get genes of interest and their gene IDs
interest_genes = np.loadtxt("dicty/Data/2017/jump_genes.csv", delimiter = ',', dtype = str, encoding="utf_8_sig")
ids = pd.read_csv("dicty/Data/updated_gene_IDs.csv", header=None)

# Remove genes beginning with DDB
interest_genes = interest_genes[~np.char.startswith(interest_genes, 'DDB')]

# Create dictionary mapping gene name to ID
gene_ids = []
gene_name_to_id = pd.Series(ids[0].values, index=ids[1]).to_dict()
for i in range(len(interest_genes)):
    gene_ids.append(gene_name_to_id[interest_genes[i]])
gene_ids = np.array(gene_ids)

# Get row indices of genes of interest
gene_indices = [np.where(rates[:, 0]== item)[0][0] for item in gene_ids]
rates = rates[gene_indices, 1].astype(float)
rates = np.insert(rates, 0, 1)
rates = np.column_stack((rates, rates/3))


np.savetxt("dicty/Rates/degradation_rates.txt", rates, delimiter = '\t', fmt = '%f')
