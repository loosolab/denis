import numpy as np
import json
import pandas as pd

# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial import distance
from tobias.utils.motifs import MotifList

def make_MotifList(file_paths):
    '''
    Create a MotifList object from list of motif files.
    '''
    # load motifs
    motif_list = MotifList()
    
    for file in file_paths:
        motif_list = motif_list + MotifList().from_file(file)
    
    return MotifList(motif_list)

##### load motifs
discovered_motifs = make_MotifList(snakemake.input.d_motifs)
control_motifs = make_MotifList(snakemake.input.c_motifs)
db_motifs = MotifList().from_file(snakemake.input.db)

if len(discovered_motifs) < 1:
    raise Exception("MotifList is empty!")

discovered_motifs_names = [f"{m.id} {m.name}" for m in discovered_motifs]
control_motifs_names = [f"{m.id} {m.name}" for m in control_motifs]
db_motifs_names = [f"{m.id} {m.name}" for m in db_motifs]

all_motifs = MotifList(discovered_motifs + control_motifs + db_motifs)

##### cluster
all_motifs.cluster(threshold=0.3, metric="pcc", clust_method="average")

# filter similarity matrix
new_vs_db = all_motifs.similarity_matrix.loc[discovered_motifs_names + control_motifs_names, db_motifs_names]

##### plot
## plot new motifs vs database heatmap
col_dendrogram_height = 2 / (len(new_vs_db) * 0.25)
plot = sns.clustermap(new_vs_db, 
                      cmap="YlOrRd_r",
                      vmin=0,
                      vmax=1,
                      xticklabels=False,
                      yticklabels=True,
                      figsize=(10, len(new_vs_db) * 0.25),
                      dendrogram_ratio=(0.2, col_dendrogram_height),
                      cbar_pos=(0.02, 1 - col_dendrogram_height, 0.05, col_dendrogram_height)
                      )
plot.fig.suptitle("Distance heatmap") 
plot.ax_heatmap.set_xlabel("Database Motifs")
plot.ax_heatmap.set_ylabel("New Motifs")
plot.savefig(snakemake.output.heatmap, bbox_inches='tight', dpi=100)
plt.close()

## plot lineplot for each new motif
db_vs_new = new_vs_db.T

n_cols = 2
n_rows = int(np.ceil(len(new_vs_db) / n_cols))

fig, axs = plt.subplots(ncols=n_cols, nrows=n_rows, figsize=(14, 2 * n_rows), tight_layout={'rect': (0, 0, 1, 0.95)})
fig.suptitle('New motifs vs. distance sorted motif database', fontsize=16)
flat_axs = axs.flatten()

top_x = dict()
for i, motif in enumerate(db_vs_new.columns):
    cur_motif_top = db_vs_new.nsmallest(snakemake.params.top_x, motif)[motif]
    top_x[motif] = (cur_motif_top.index.tolist(), cur_motif_top.values.tolist())
    
    ##### plot #####
    plot = sns.lineplot(y=sorted(db_vs_new[motif]), x=range(len(db_vs_new)), ax=flat_axs[i])
    plot.set(ylim=(0, 1), ylabel='Distance')
    plot.set_title(motif)
    
    # place a text box in upper left in axes coords
    txt = "\n".join([f"Min value: {np.min(db_vs_new[motif])}", f"Max value: {np.max(db_vs_new[motif])}"])
    plot.text(0.8, 0.08, txt, transform=plot.transAxes, fontsize=8, verticalalignment='bottom')

# delete empty subplots
if i < len(flat_axs):
    for j in range(i+1, len(flat_axs)):
        fig.delaxes(flat_axs[j])

fig.savefig(snakemake.output.lineplot, dpi=100)
plt.close()

## plot distance boxplot
plt.figure(figsize=(len(db_vs_new.columns) * 0.25, 10))
boxplot = sns.boxplot(data=db_vs_new)

boxplot.set(ylabel='Distance')
boxplot.set_title('Distance Overview')
boxplot.set_xticklabels(boxplot.get_xticklabels(),rotation=90)

boxplot.get_figure().savefig(snakemake.output.boxplot, bbox_inches='tight', dpi=100)
plt.close()

# save topx most similar motifs for each new motif
with open(snakemake.output.json, "w") as file:
    json.dump(top_x, file, indent=4)
