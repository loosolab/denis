# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# load naming table
name_table = pd.read_csv(snakemake.input.motif_names, sep="\t")

# construct feature table
feature_table = pd.DataFrame(columns=["motif", "feature", "enrichment"])

for f_table in snakemake.input.hits_tables:
    hits = pd.read_csv(f_table, sep="\t")
    
    filename = hits.loc[0, "peak_id"]
    motif_name = name_table.loc[name_table["filafy_id_name"] == filename, "id_name"].values[0]
    
    bed = pd.read_csv(os.path.join(snakemake.input.bed_tables_path, f"{filename}.bed"), sep="\t", header=None)
    
    # compute enrichment
    # compute unique_peaks as a peak can be annotated multiple times
    unique_peaks = len(bed)
    counts = hits["feature"].value_counts(dropna=False)
    enrichment = counts / unique_peaks
    
    motif_table = pd.DataFrame({"motif": motif_name, "feature": enrichment.index, "enrichment": enrichment.values, "nsites_bound": unique_peaks})
    
    feature_table = feature_table.append(motif_table, ignore_index=True, sort=False)

# set feature NaN to 'Other'
feature_table.loc[feature_table["feature"].isnull(), "feature"] = "Other"

# sort
feature_table.sort_values(by=["feature", "enrichment"], inplace=True)

# add enrichment to ranks table
ranks_table = pd.read_csv(snakemake.input.ranks_table, sep="\t")
ranks_table["tmp_id"] = ranks_table["id"] + "_" + ranks_table["name"]

# rearrange to wide format
wide_features = feature_table.pivot(columns="feature", index=["motif", "nsites_bound"], values="enrichment")
# squash multilevel index
wide_features.reset_index(inplace=True)

# set NA features to 0
wide_features.fillna(0, inplace=True)

# merge
ranks_table = ranks_table.merge(wide_features, how="left", left_on="tmp_id", right_on="motif")
ranks_table.drop(["tmp_id", "motif"], axis=1, inplace=True)

# save table
ranks_table.to_csv(snakemake.output.table, sep="\t", index=False)

# melt table (wide to long format) for plotting
feature_table = ranks_table.melt(id_vars=["id", "name", "GC%", "nsites_whole_genome", "nsites_open_chromatin", "nsites_bound", "width", "information_content", "enrichment_score"], var_name="feature", value_name="feature_enrichment")

# plot
plt.figure(figsize=(len(snakemake.input.hits_tables) * 0.25, 10))

# use seaborn theme (enable gridlines)
sns.set_theme()

xaxis = "name"
plot = sns.scatterplot(data=feature_table,
                       x=xaxis,
                       y="feature_enrichment",
                       hue="feature",
                       size="nsites_bound"
)

# https://stackoverflow.com/questions/66700765/how-to-rotate-x-axis-tick-lables-in-seaborn-scatterplot-using-subplots
plot.set_xticklabels(feature_table[xaxis].to_list(), rotation=90)
# force x-gridlines
plot.xaxis.grid(True)

plot.get_figure().savefig(snakemake.output.plot, bbox_inches='tight', dpi=100)
plt.close()
