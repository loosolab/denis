import os
import pybedtools
import pandas as pd

# load all extracted footprints
footprint_bed = pybedtools.BedTool(snakemake.input.footprint_bed)
footprint_bed = footprint_bed.sort()

# overlap motif rescan bed with all footprints
# load bed file
rescan = pybedtools.BedTool(snakemake.input.rescan_bed)
rescan = rescan.sort()

# keep all footprints that overlap
overlap = footprint_bed.intersect(rescan, wa=True)

# add name
name = os.path.basename(os.path.splitext(snakemake.input.rescan_bed)[0])
overlap_table = overlap.to_dataframe()
overlap_table["name"] = name

# save
overlap_table.to_csv(snakemake.output.bed, sep="\t", index=False, header=None)
