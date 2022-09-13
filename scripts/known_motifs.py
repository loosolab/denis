# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import os
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial import distance
from tobias.utils.motifs import MotifList

def get_cluster_reps(linkage, t, label):
    '''
    Get the middle element of each cluster in a dendrogram.
    
    Parameters:
    linkage (numpy.array): Output of scipy.cluster.hierarchy.linkage
    t (int): Number of returned representatives (See scipy.cluster.hierarchy.fcluster)
    label (list): List of labels according to original data used to create linkage.
    
    Returns:
    List of cluster representatives
    '''
    # compute clusters
    clusters = hierarchy.fcluster(linkage, t, criterion='maxclust')
    
    dn = hierarchy.dendrogram(linkage, labels=label)
    plt.close()
    
    # get cluster members
    cluster_dict = {c: [label[i] for i, c2 in enumerate(clusters) if c == c2] for c in set(clusters)}
    # order members like dendrogram
    cluster_dict = {key: [e for e in dn["ivl"] if e in value] for key, value in cluster_dict.items()}

    # select middle motif
    reps = [values[int((len(values) - 1) / 2)] for values in cluster_dict.values()]

    return reps

# create output dir if necessary
os.makedirs(snakemake.output[0], exist_ok=True)

# skip everything if no motif file was provided
if snakemake.params.file:
    # load motifs
    motifs = MotifList().from_file(snakemake.params.file)

    motif_names = []
    ## select motifs based on distance
    if snakemake.params.quantity:
        # cluster motifs based on similarity
        motifs.cluster(threshold=0.3, metric="pcc", clust_method="average")

        # select representatives
        linkage = hierarchy.linkage(distance.pdist(motifs.similarity_matrix), method="average")
        motif_names += get_cluster_reps(linkage, snakemake.params.quantity, motifs.similarity_matrix.index)

    ## select motifs from custom selection in config
    if snakemake.params.motifs:
        all_names = [(m.id, m.name) for m in motifs]

        motif_names += [f"{id} {name}" for id, name in all_names if id in snakemake.params.motifs or name in snakemake.params.motifs]

    # save motifs and logos
    motif_obj = [m for m in motifs if f"{m.id} {m.name}" in motif_names]

    for m in motif_obj:
        m.to_file(os.path.join(snakemake.output[0], f"{m.id}_{m.name}.meme"), fmt="meme")
        m.logo_to_file(os.path.join(snakemake.output[0], f"{m.id}_{m.name}.png"))
