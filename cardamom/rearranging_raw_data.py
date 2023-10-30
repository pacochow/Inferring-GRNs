import numpy as np

# Get expression matrix for all genes
full_data = np.loadtxt("dicty/Data/2023/raw_data.csv", delimiter = ',', dtype = str)

norm_data_split = np.loadtxt("dicty/Data/2023/divided_with_labels.csv", delimiter = ',', dtype = str).astype('<U14')
def rearrange_columns_to_match(A, B):
    # Use the first row as a representative to find the column order
    representative_row_A = np.array([item.strip('"') for item in A[0]])
    representative_row_B = B[0]

    # Get the indices that would sort A's row to match B's row
    order = [np.where(representative_row_A == item)[0][0] for item in representative_row_B]

    # Re-arrange the columns of A using the order obtained
    A_reordered = A[:, order]

    return A_reordered

full_data = rearrange_columns_to_match(full_data, norm_data_split)
np.savetxt("dicty/Data/2023/rearranged_raw_data.csv", full_data, delimiter = ',', fmt = '%s')
