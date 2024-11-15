import os
from datetime import datetime

start_time = datetime.now().strftime("%a %b %-d %H:%M:%S %Y")

config: 'config.yml'

#-------------------------------------------------------------------------------#
#--------------------------- Generate output files -----------------------------#
#-------------------------------------------------------------------------------#

OUTPUTDIR = config['run_params']['output']

output_files = ["config.yml",
                "2_discovery/2_processed_motif/done.txt",
                "2_discovery/motif_stats.tsv",
                "3_evaluation/motif_ranks.tsv",
                "3_evaluation/ranks_barplot.pdf"]

# skip some evaluation steps if no motif database is provided
if config["footprint"]["motif_file"]:
    output_files.extend(
        ["3_evaluation/motif_evaluation/new_vs_db_heatmap.pdf",
         "3_evaluation/motif_evaluation/new_vs_db_lineplot.pdf",
         "3_evaluation/motif_evaluation/new_vs_db_boxplot.pdf",
         "3_evaluation/motif_evaluation/most_similar_motifs.json"]
    )

if config['annotation']['gtf']:
    # annotation files
    output_files.append("4_annotation/feature_enrichment_plot.pdf")
    output_files.append("4_annotation/feature_enrichment_table.tsv")

    if config['go_enrichment']['email']:
        # gene ontology files
        output_files.append("5_gene_ontology/done.txt")

output_files = [os.path.join(OUTPUTDIR, f) for f in output_files]

#-------------------------------------------------------------------------------#
#---------------------------------- RUN :-) ------------------------------------#
#-------------------------------------------------------------------------------#

include: 'snakefiles/helpers.snake'
include: 'snakefiles/footprint.snake'
include: 'snakefiles/motif.snake'
include: 'snakefiles/evaluation.snake'
include: 'snakefiles/annotation.snake'
include: 'snakefiles/gene_ontology.snake'

rule all:
    input:
        output_files
    message: f"Pipeline done. Started {start_time}"
