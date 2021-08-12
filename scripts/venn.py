# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import pybedtools
import matplotlib_venn as mv
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import numpy as np
import seaborn as sns
import pandas as pd
import warnings

# create tmp folder for pybedtools
tmp_dir = os.path.join(snakemake.output.dir, "tmp")
os.makedirs(tmp_dir, exist_ok=True)
pybedtools.helpers.set_tempdir(tmp_dir)

def venn(set_, group_names, legend_title, ax, scientific_notation=True, subset_label=True, legend_label=None):
    # https://github.com/konstantint/matplotlib-venn/issues/30
    # catch layout error by using equal proportions
    try:
        plot = mv.venn3(subsets = set_, set_labels=group_names if subset_label else None, ax=ax)
    except mv._region.VennRegionException:
        # equal_size = {'001': 1, '010': 1, '011': 1, '100': 1, '101': 1, '110': 1, '111': 1}
        equal_size = {key: 100 if value > 0 else value for key, value in set_.items()}
        plot = mv.venn3(subsets = equal_size, set_labels=group_names if subset_label else None, ax=ax)

    colors = {'001': '#9999FF', 
              '010': '#99CD99', 
              '011': 'lightblue', 
              '100': '#FF9999', 
              '101': '#E199E1', 
              '110': '#E1BD99', 
              '111': '#C2AEC2'}
    # https://stackoverflow.com/a/46224448
    h, l = [],[]
    for i, key in enumerate(set_):
        # remove label by setting them to empty string:
        if not plot.get_label_by_id(key) is None:
            plot.get_label_by_id(key).set_text("")
        
        patch = mpatches.Patch() if plot.get_patch_by_id(key) is None else plot.get_patch_by_id(key)
        patch.set_color(colors[key])
        patch.set_alpha(1)
        
        # append patch to handles list
        h.append(patch)
        
        # append count to labels list
        overlap_label = legend_label[i] if legend_label else "-".join([group_names[index] for index, char in enumerate(list(key)) if char == "1"])
        if scientific_notation:
            l.append(f"{set_[key]:e} ({overlap_label})")
        else:
            l.append(f"{round(set_[key], 5)} ({overlap_label})")

    #create legend from handles and labels    
    ax.legend(handles=h, labels=l, title=legend_title, loc='lower left', bbox_to_anchor=(1, 0.5))
    
    return ax

# load peaks
peaks = pybedtools.BedTool(snakemake.input.peak_bed)
peaks = peaks.slop(g=snakemake.input.genome_file, b=snakemake.params.extend_region)
peaks = peaks.sort()

# load and prepare motif name table
motif_names = pd.read_csv(snakemake.input.motif_names, sep="\t")

# motif overview table
overview = {}

for rescan_bed, motif_bed in zip(snakemake.input.rescan_beds, snakemake.input.motif_beds):
    if not os.path.basename(rescan_bed) == os.path.basename(motif_bed):
        warnings.warn(f"Filenames {os.path.basename(rescan_bed)} and {os.path.basename(motif_bed)} don't match.", Warning)
    
    # load bed files
    rescan = pybedtools.BedTool(rescan_bed)
    footprints = pybedtools.BedTool(motif_bed)

    filename = os.path.basename(os.path.splitext(rescan_bed)[0])
    name = motif_names.loc[motif_names["filafy_id_name"] == filename, ["id", "name"]].values[0]
    id = name[0]
    full_name = name[0] + " " + name[1]
    
    # sort beds
    rescan = rescan.sort()
    footprints = footprints.sort()

    ###### prepare data #####
    # rescan
    rescan_overlap = {'001': rescan.intersect(peaks, v=True).intersect(footprints, v=True),
                       '011': rescan.intersect(peaks, wa=True, u=True).intersect(footprints, v=True),
                       '101': rescan.intersect(footprints, wa=True, u=True).intersect(peaks, v=True),
                       '111': rescan.intersect(footprints, wa=True, u=True).intersect(peaks, wa=True, u=True)}
    
    # footprint
    footprint_overlap = {'100': footprints.intersect(rescan, v=True).intersect(peaks, v=True),
                         '101': footprints.intersect(rescan, wa=True, u=True).intersect(peaks, v=True),
                         '110': footprints.intersect(peaks, wa=True, u=True).intersect(rescan, v=True),
                         '111': footprints.intersect(rescan, wa=True, u=True).intersect(peaks, wa=True, u=True)}
    
    # peak
    peak_overlap = {'010': peaks.intersect(rescan, v=True).intersect(footprints, v=True),
                    '011': peaks.intersect(rescan, wa=True, u=True).intersect(footprints, v=True),
                    '110': peaks.intersect(footprints, wa=True, u=True).intersect(rescan, v=True),
                    '111': peaks.intersect(rescan, wa=True, u=True).intersect(footprints, wa=True, u=True)}
    
    # calculate genome length
    genome_len = sum(pd.read_csv(snakemake.input.genome_file, sep="\t", header=None)[2])
    
    ##### plotting #####
    fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(17, 10), tight_layout={'rect': (0, 0, 1, 0.95)})
    axs = axs.flatten()
    fig.suptitle(f"{full_name} binding site scan", fontsize=16)
    
    ## venn ##
    plot_name = ["Rescan", "Original Sites", "Open Chromatin"]
    tick_names = ["Original Sites", "Open Chromatin", "Rescan"]
    for i, data in enumerate([rescan_overlap, footprint_overlap, peak_overlap]):
        # skip original site venn if there are no original sites
        if i == 1 and rescan_bed == motif_bed:
            fig.delaxes(axs[i])

            continue
        
        # count overlaps and set legend labels
        count_data = {}
        legend_label = []
        for key, value in data.items():
            count_data[key] = value.count()
            
            if key.count("1") == 1:
                legend_label.append("No Overlap")
            else:
                legend_label.append(" and ".join([tick_names[index] for index, char in enumerate(list(key)) if char == "1" and tick_names[index] != plot_name[i]]))
        
        try:
            venn(count_data, group_names=tick_names, legend_title=f"{plot_name[i]} Overlap Count", ax=axs[i], scientific_notation=False, subset_label=False, legend_label=legend_label)
        except mv._region.VennRegionException:
            # make alternative bar plot
            counts = pd.DataFrame({"Label": ["-".join([tick_names[index] for index, char in enumerate(list(k)) if char == "1"]) for k in count_data.keys()],
                                "Count": count_data.values()})

            sns.barplot(x="Label", y=f"{plot_name[i]} Overlap Count", data=counts, ax=axs[i])

            axs[i].set_yscale('log')
            axs[i].get_yaxis().set_major_formatter(ScalarFormatter())
            axs[i].ticklabel_format(style='plain', axis='y')
            
            axs[i].set_xlabel('')
            axs[i].set_xticklabels(axs[i].get_xticklabels(), rotation=90)

    ## barplot ##
    # enrichment scores
    enrichment = pd.DataFrame({"Label": ["Genome Enrichment", "Open Chromatin Enrichment"], 
                               "Enrichment": [rescan.total_coverage() / genome_len * 100 ,
                               (rescan_overlap['011'].total_coverage() + rescan_overlap['111'].total_coverage()) / peaks.total_coverage() * 100]})
    
    sns.barplot(x="Label", y="Enrichment", data=enrichment, ax=axs[3])
    axs[3].set_xlabel('')
    axs[3].set_ylabel('Basepair Ratio [%]')
    
    plt.savefig(os.path.join(snakemake.output.dir, filename + ".pdf"))
    plt.close()
    
    # add to overview table
    overview.setdefault("id", []).append(id)
    overview.setdefault("enrichment_score", []).append(np.log2(enrichment.loc[1, "Enrichment"] / enrichment.loc[0, "Enrichment"]))
    overview.setdefault("nsites_whole_genome", []).append(len(rescan))
    overview.setdefault("nsites_open_chromatin", []).append(rescan.intersect(peaks, wa=True, u=True).count())

# convert to table and write
overview = pd.DataFrame(overview).sort_values(by="enrichment_score", ascending=False)
overview.to_csv(snakemake.output.ranks, sep="\t", index=False)
