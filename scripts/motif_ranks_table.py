# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import pandas as pd
from tobias.utils.motifs import MotifList
import seaborn as sns
import matplotlib.pyplot as plt

# load motifs
motifs = MotifList().from_file(snakemake.input.motifs)

# create motif info dict
motif_dict = {}

for m in motifs:
    m.information_content()
    m.gc_content()
    
    motif_dict.setdefault("id", []).append(m.id)
    motif_dict.setdefault("name", []).append(m.name)
    motif_dict.setdefault("GC%", []).append(m.gc)
    motif_dict.setdefault("nsites", []).append(int(m.info['nsites']))
    motif_dict.setdefault("width", []).append(int(m.info['w']))
    motif_dict.setdefault("information_content", []).append(m.ic)

motif_table = pd.DataFrame(motif_dict)
motif_table["merge_col"] = motif_table["id"] + "_" + motif_table["name"]

# load rank table
ranks = pd.read_csv(snakemake.input.table, sep="\t")
ranks.rename(columns={"id": "merge_col"}, inplace=True)

# merge tables
ranks = motif_table.merge(ranks, on="merge_col", how="left")
ranks.drop("merge_col", axis=1, inplace=True)
ranks.sort_values("enrichment_score", axis=0, ascending=False, inplace=True)

# write table
ranks.to_csv(snakemake.output[0], sep="\t", index=False)

# rank barplot
plt.figure(figsize=(len(ranks) * 0.25, 10))

plot = sns.barplot(data=ranks, x=ranks["id"] + " " + ranks["name"], y="enrichment_score")
plot.set_xticklabels(plot.get_xticklabels(), rotation=90)
plot.set(xlabel='', ylabel='Enrichment Score')

#save plot
plot.get_figure().savefig(snakemake.output[1], bbox_inches='tight', dpi=100)
plt.close()
