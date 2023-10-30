import numpy as np
from tqdm import tqdm

# Get expression matrix for all genes
categorised_cells = np.loadtxt("dicty/Data/2023/categorised_cells.csv", delimiter = ',', dtype = str)


# Get expression matrix for all genes
full_data = np.loadtxt("dicty/Data/2023/Norm_data_SCRAN.csv", delimiter = ',', dtype = str)
full_data[0] = np.array([i.strip('"') for i in full_data[0]])

processed_data = np.zeros((len(full_data), len(categorised_cells)+1)).astype(str)
processed_data[:, 0] = full_data[:, 0]
processed_data[0, 1:] = categorised_cells[:, 1]
for i in tqdm(range(len(categorised_cells))):
    cell_name = categorised_cells[i,0]
    full_data_index = np.where(full_data[0] == cell_name)[0][0]
    processed_data[1:, i+1] = full_data[1:, full_data_index]
    
np.savetxt("dicty/Data/2023/clustered_data_norm.csv", processed_data, delimiter = ',', fmt = '%s')