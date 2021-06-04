# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# construct feature table
feature_table = pd.DataFrame(columns=["motif", "feature", "enrichment"])

for f_table in snakemake.input.hits_tables:
    hits = pd.read_csv(f_table, sep="\t")
    
    motif_name = hits.loc[0, "peak_id"]
    
    # compute enrichment
    # compute unique_peaks as a peak can be annotated multiple times
    unique_peaks = len(set(hits.iloc[:,0] + hits.iloc[:,1].astype(str) + hits.iloc[:,2].astype(str)))
    counts = hits["feature"].value_counts(dropna=False)
    enrichment = counts / unique_peaks
    
    motif_table = pd.DataFrame({"motif": motif_name, "feature": enrichment.index, "enrichment": enrichment.values})
    
    feature_table = feature_table.append(motif_table, ignore_index=True, sort=False)

# set feature NaN to 'Other'
feature_table.loc[feature_table["feature"].isnull(), "feature"] = "Other"

# sort
feature_table.sort_values(by=["feature", "enrichment"], inplace=True)

# plot
plt.figure(figsize=(len(snakemake.input.hits_tables) * 0.25, 10))

# use seaborn theme (enable gridlines)
sns.set_theme()
plot = sns.pointplot(data=feature_table,
                     x="motif",
                     y="enrichment",
                     hue="feature")

plot.set_xticklabels(plot.get_xticklabels(), rotation=90)
# force x-gridlines
plot.xaxis.grid(True)

plot.get_figure().savefig(snakemake.output.plot, bbox_inches='tight', dpi=100)
plt.close()

# add enrichment to ranks table
ranks_table = pd.read_csv(snakemake.input.ranks_table, sep="\t")
ranks_table["tmp_id"] = ranks_table["id"] + "_" + ranks_table["name"]

# rearrange to wide format
wide_features = feature_table.pivot(columns="feature", index="motif", values="enrichment")

# merge
ranks_table = ranks_table.merge(wide_features, how="left", left_on="tmp_id", right_on="motif")
ranks_table.drop("tmp_id", axis=1, inplace=True)

# save table
ranks_table.to_csv(snakemake.output.table, sep="\t", index=False)
