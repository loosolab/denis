from Bio import SeqIO

with open(snakemake.output.genome_file, "w") as output:
    for record in SeqIO.parse(snakemake.input.genome, "fasta"):
        # add custom third column; length without 'N' bases
        output.write(f"{record.id}\t{len(record)}\t{len(record) - record.upper().seq.count('N')}\n")
