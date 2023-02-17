# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import os
from pathlib import Path
from tobias.utils.motifs import MotifList
from tobias.tools.motif_clust import plot_dendrogram
import pandas as pd

# create motif list
motif_list = MotifList()

for meme_file in snakemake.input.motifs:
    motif_list += MotifList().from_file(meme_file)
    
    # save source file
    motif_list[-1].source_file = meme_file

if len(motif_list) < 1:
    raise Exception("MotifList is empty!")

# cluster then save distance matrix
cluster = motif_list.cluster(threshold=snakemake.params.threshold, metric=snakemake.params.metric, clust_method="average")

# create path if needed
motif_list.similarity_matrix.to_csv(os.path.join(snakemake.output.dir, "distance.tsv"), sep = '\t')
# plot dendrogram
plot_dendrogram(motif_list.similarity_matrix.columns, 
                motif_list.linkage_mat,
                font_size=12,
                out=os.path.join(snakemake.output.dir, "dendrogram.png"),
                title="Clustering",
                threshold=snakemake.params.threshold,
                dpi=snakemake.params.dpi)

# output path for consensus motifs
out_path = os.path.join(snakemake.output.dir, "motifs")
os.makedirs(out_path, exist_ok=True)

# load motif stats table
motif_stats = pd.read_csv(snakemake.input.stats, sep="\t")
motif_stats["consensus"] = None

# create consensus motifs by combining motifs of a cluster
# combine motifs of each cluster to consensus
consensus_motifs = MotifList()
for i, cluster_id in enumerate(cluster):
    consensus = cluster[cluster_id].create_consensus(metric=snakemake.params.metric) # MotifList object with create_consensus method
    
    consensus_name = f"motif_{i}"
    consensus.id = consensus_name # cluster_id if len(cluster[cluster_id]) > 1 else cluster[cluster_id][0].id # set original motif id if cluster length = 1
    consensus.name = consensus_name
    
    # TODO fix add length because of bug
    consensus.length = len(consensus.counts[0])
    
    site_count = 0
    # save consensus motif sites
    # https://stackoverflow.com/a/17749339
    with open(os.path.join(out_path, f"{consensus.id}_{consensus.name}.fasta"), "w") as out_fasta:
        with open(os.path.join(out_path, f"{consensus.id}_{consensus.name}.bed"), "w") as out_bed:
            for motif in cluster[cluster_id]:
                # add name to motif stats
                motif_stats.loc[(motif_stats["id"] == motif.id) & (motif_stats["name"] == motif.name), "consensus"] = consensus_name
                
                # replace file extension with .fasta
                input_fasta = os.path.splitext(motif.source_file)[0] + ".fasta"
                
                with open(input_fasta, "r") as in_fasta:
                    for line in in_fasta:
                        # write fasta
                        out_fasta.write(line)
                        
                        # write bed
                        if line.startswith(">"):
                            _, _, chr_, start, end = line.rstrip().replace("-", ":").split(":")
                            out_bed.write("\t".join([chr_, start, end, consensus.id]) + "\n")
                            
                            site_count +=1
    
    # add site count to motif
    consensus.n = site_count

    consensus_motifs.append(consensus)

# Save motif stats
motif_stats.to_csv(snakemake.output.stats, sep="\t", index=False)

# Save meme and logo of every motif
for motif in consensus_motifs:
    # save meme file
    motif.to_file(os.path.join(out_path, f"{motif.id}_{motif.name}.meme"), "meme")
    
    # save logo png
    filename = os.path.join(out_path, f"{motif.id}_{motif.name}.png")
    motif.logo_to_file(filename)

Path(snakemake.output.done).touch()
