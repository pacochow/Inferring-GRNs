from src.analysis_utils import *
import markov_clustering as mc

for i in [None, 0, 1, 2]:
    draw_graph(2017, i, 2, False)

# inter = np.load("dicty/cardamom/inter.npy")
# G, pos = create_graph(inter)
# # Run MCL with default parameters
# result = mc.run_mcl(inter)           # Run MCL with default parameters
# clusters = mc.get_clusters(result)   # Get clusters

# # Print clusters
# for cluster in clusters:
#     print(cluster)


# new_inter = extract_top_interactions(100, inter)
# analytics = compute_analytics(new_inter)
# analytics = analytics.sort_values(by = 'Weighted score')
# analytics['Rank'] = analytics['Weighted score'].rank(ascending = False).astype(int)
# print(analytics.sort_values(by = 'Rank'))

# # Spectral Clustering
# n_clusters = 5  # Specify the number of clusters you want
# sc = SpectralClustering(n_clusters=n_clusters, affinity='precomputed')
# labels = sc.fit_predict(inter)

