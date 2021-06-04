from tobias.utils.motifs import MotifList

motif_list = MotifList()

for file in snakemake.input:
    motif_list += MotifList().from_file(file)

motif_list.to_file(snakemake.output[0], "meme")
