import pandas as pd

def gene_id_dicts():
    
    # Get array of gene names corresponding to gene IDs
    ids = pd.read_csv("dicty/Data/updated_gene_IDs.csv", header=None)

    # Create dictionary mapping gene names to gene IDs
    gene_ids = []
    gene_name_to_id = pd.Series(ids[0].values, index=ids[1]).to_dict()
    gene_id_to_name = pd.Series(ids[1].values, index = ids[0]).to_dict()
    
    return gene_name_to_id, gene_id_to_name


def jump_gene_category_dicts():
    
    jump_genes = pd.read_csv("dicty/Data/jump_genes_categorised.csv", header = None)
    
    gene_name_to_cat = pd.Series(jump_genes[2].values, index=jump_genes[1]).to_dict()
    
    return gene_name_to_cat