import os
from datetime import datetime

start_time = datetime.now().strftime("%a %b %-d %H:%M:%S %Y")

config: 'config.yml'

#-------------------------------------------------------------------------------#
#--------------------------- Generate output files -----------------------------#
#-------------------------------------------------------------------------------#

OUTPUTDIR = config['run_params']['output']

output_files = ["2_discovery/2_processed_motif/done.txt",
                "2_discovery/motif_stats.tsv",
                "3_evaluation/motif_evaluation/new_vs_db_heatmap.pdf",
                "3_evaluation/motif_evaluation/new_vs_db_lineplot.pdf",
                "3_evaluation/motif_evaluation/new_vs_db_boxplot.pdf",
                "3_evaluation/motif_evaluation/most_similar_motifs.json",
                "3_evaluation/motif_ranks.tsv",
                "3_evaluation/ranks_barplot.pdf",
                "config.yml"]

if config['annotation']['gtf']:
    output_files.append("4_annotation/feature_enrichment_plot.pdf")
    output_files.append("4_annotation/feature_enrichment_table.tsv")

output_files = [os.path.join(OUTPUTDIR, f) for f in output_files]

#-------------------------------------------------------------------------------#
#---------------------------------- RUN :-) ------------------------------------#
#-------------------------------------------------------------------------------#

include: 'snakefiles/footprint.snake'
include: 'snakefiles/motif.snake'
include: 'snakefiles/evaluation.snake'
include: 'snakefiles/annotation.snake'
include: 'snakefiles/helpers.snake'

rule all:
    input:
        output_files
    message: f"Pipeline done. Started {start_time}"