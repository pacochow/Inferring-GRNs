import numpy as np
import pandas as pd
import networkx as nx
import copy
import matplotlib.pyplot as plt
from pyvis.network import Network
from src.processing_utils import *

def get_used_genes():
    return np.loadtxt("dicty/Data/panel_genes.txt", delimiter = '\t', dtype = str, encoding="utf_8_sig")

def get_indices_names_dicts(genes):
    indices_to_gene_names = {int(row[0]): row[1] for row in genes}
    gene_name_to_indices = {row[1]: int(row[0]) for row in genes}
    return indices_to_gene_names, gene_name_to_indices
    
def get_scores(inter, gene_index):
    # Get interaction matrix

    outgoing = (inter[gene_index] != 0).sum()
    incoming = (inter[:, gene_index] != 0).sum()
    score = np.abs(inter[gene_index]).sum()+np.abs(inter[:,gene_index]).sum()-np.abs(inter[gene_index,gene_index])
    
    return outgoing, incoming, score
    
def extract_top_interactions(percent_cutoff, inter):
    new_inter = copy.deepcopy(inter)
    n_positive_interactions = (inter>0).sum()
    n_negative_interactions = (inter<0).sum()
    # Extract top % interactions
    for i in range(len(inter)):
        for j in range(len(inter)):
            if inter[i,j] > 0 and inter[i,j]<np.sort(inter, axis=None)[::-1][int(percent_cutoff*0.01*n_positive_interactions)]:
                new_inter[i,j] = 0
            elif inter[i,j] < 0 and inter[i,j] > np.sort(inter, axis=None)[int(percent_cutoff*0.01*n_negative_interactions)]:
                new_inter[i,j] = 0
    return new_inter
    
def compute_analytics(inter):
    indices_to_gene_name, _ = get_indices_names_dicts(get_used_genes())
    
    analytics = pd.DataFrame(columns = ('Gene name', 'Outgoing connections', 'Incoming connections', 'Weighted score'))
    for gene_index in range(len(inter)):
        outgoing, incoming, score = get_scores(inter, gene_index)
        analytics.loc[gene_index] = [indices_to_gene_name[gene_index], outgoing, incoming, score]    
    return analytics

def create_graph(inter):
    indices_to_gene_name, _ = get_indices_names_dicts(get_used_genes())
    # Convert the numpy matrix to a directed graph
    G = nx.from_numpy_array(inter, create_using=nx.DiGraph)
    pos_original = nx.shell_layout(G)

    # Rename nodes
    G_renamed = nx.relabel_nodes(G, indices_to_gene_name)

    # Update the positions dictionary to use the new node names
    pos_renamed = {indices_to_gene_name[old_name]: pos for old_name, pos in pos_original.items()}
    return G_renamed, pos_renamed

    
def analyse_graph(inferred_timepoint_index):

    
    if inferred_timepoint_index == None:
        inter = np.load("dicty/cardamom/inter.npy")
    else:
        inter_t = np.load("dicty/cardamom/inter_t.npy")
        inter = inter_t[inferred_timepoint_index]
        
    for cutoff in np.arange(1, 10):
        new_inter = extract_top_interactions(cutoff, inter)
       
def draw_graph(dataset, inferred_timepoint_index = None, percent_cutoff = 2, interactive_graph = True):

    # Get timepoint labels
    if dataset == 2023:
        labels = [0, 1, 7, 9]
    elif dataset == 2017:
        labels = [0, 3, 6]

    # Get genes of interest and their gene IDs
    interest_genes = np.loadtxt(f"dicty/Data/{dataset}/jump_genes.csv", delimiter = ',', dtype = str, encoding="utf_8_sig")
    ids = pd.read_csv("dicty/Data/updated_gene_IDs.csv", header=None)

    # Remove genes beginning with DDB
    interest_genes = interest_genes[~np.char.startswith(interest_genes, 'DDB')]

    # Create dictionary mapping gene IDs to gene names
    _, gene_id_to_name = gene_id_dicts()
    gene_id_to_name["Stimulus"] = "Stimulus"

    # Create dictionary mapping gene IDs to categories
    gene_name_to_cat = jump_gene_category_dicts()
    gene_name_to_cat["Stimulus"] = "Stimulus"

    # Get interaction matrix
    if inferred_timepoint_index == None:
        inter = np.load("dicty/cardamom/inter.npy")
        label = "final"
    else:
        inter_t = np.load("dicty/cardamom/inter_t.npy")
        inter = inter_t[inferred_timepoint_index]
        label = labels[inferred_timepoint_index]

    # Create indices to gene name dictionary
    genes = get_used_genes() 
    indices_to_gene_name = {int(row[0]): row[1] for row in genes}

    colour_map = ['orange' if gene_name_to_cat[indices_to_gene_name[i]] == "Pre-jump" else 'skyblue' if gene_name_to_cat[indices_to_gene_name[i]] == "Jump" 
                  else 'violet' if gene_name_to_cat[indices_to_gene_name[i]] == "Post-jump" else 'Gray' for i in indices_to_gene_name]


    inter = extract_top_interactions(percent_cutoff, inter)

    G_renamed, pos_renamed = create_graph(inter)

    

    # Increase the figure size
    plt.figure(figsize=(11, 11))

    # Draw nodes and labels
    nx.draw_networkx_nodes(G_renamed, pos_renamed, node_size=800, node_color=colour_map, alpha=0.6)
    nx.draw_networkx_labels(G_renamed, pos_renamed, font_size=8)

    # Extract edges with positive and negative weights
    positive_edges = [(u, v) for u, v, d in G_renamed.edges(data=True) if d['weight'] > 0]
    negative_edges = [(u, v) for u, v, d in G_renamed.edges(data=True) if d['weight'] < 0]

    weights = [G_renamed[u][v]['weight'] for u,v in G_renamed.edges()]
    positive_weights = (np.array(weights)[[i>0 for i in weights]])
    negative_weights = (np.array(weights)[[i<0 for i in weights]])

    # Draw positive edges with blue color and negative edges with red color
    nx.draw_networkx_edges(G_renamed, pos_renamed, edgelist=positive_edges, edge_color="green", width=10*(positive_weights/np.sum(weights)))
    nx.draw_networkx_edges(G_renamed, pos_renamed, edgelist=negative_edges, edge_color="red", width=10*(negative_weights/np.sum(weights)))

    plt.savefig(f"dicty/Results/{dataset}_norm_graph_{label}.png", bbox_inches = 'tight')
    plt.show()

    if interactive_graph == True:
        # Create an interactive plot
        nt = Network(notebook=True, width="100%", height="800px")
        # Add nodes
        for i, node in enumerate(G_renamed.nodes()):
            nt.add_node(node, color = colour_map[i])

        # Add edges with color, thickness, and title based on weight
        for u, v, data in G_renamed.edges(data=True):
            color = 'green' if data['weight'] > 0 else 'red'
            title = f"Weight: {data['weight']:.3f}"  # Display weight as title
            nt.add_edge(u, v, color=color, value=np.exp(data['weight']**2), title=title)


        # Emphasize directed edges
        for edge in nt.edges:
            edge['arrows'] = 'to'
            edge['arrowStrikethrough'] = True  # This makes the arrow appear on top of the edge
            edge['scaling'] = {
                'min': 1,
                'max': 3,
                'label': True
            }

        # Adjust physics settings for barnesHut solver
        nt.physics = {
            "enabled": True,
            "solver": "barnesHut",
            "barnesHut": {
                "gravitationalConstant": -1,  # Controls the repulsion between nodes
                "centralGravity": 0.1,           # Controls the attraction towards the center
                "springLength": 500,             # Natural length of the springs (edges)
                "springConstant": 0.01,          # Controls the springiness of the edges
                "damping": 1,                  # Higher damping leads to more stable layout
                "avoidOverlap": 0.1              # Tries to avoid node overlap
            }
        }

        nt.show(f"dicty/Results/{dataset}_norm_graph_{label}.html")