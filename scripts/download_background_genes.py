import pandas as pd
from Bio import Entrez
import time

def chunks(lst, n):
    """
    Yield successive n-sized chunks from lst.
    https://stackoverflow.com/a/312464/19870975
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# ----- identify organism taxid ----- #
# load gtf
gtf = pd.read_csv(snakemake.input.gtf, sep="\t", header=None)

# get gene id to identify organism
for gene in gtf[8].str.extract(r'gene_id "(.+?)".*$')[0]:
    if str(gene).startswith("ENS"):
        break

# load gene2ensembl
# contains taxid and entrez to ensembl id mapping
gene2ensembl = pd.read_csv(snakemake.input.gene2ensembl, sep="\t")

taxid = gene2ensembl.loc[gene2ensembl["Ensembl_gene_identifier"] == gene, "#tax_id"].values[0]

# ----- get protein coding genes ----- #
# ncbi search term
term = f'"{taxid}"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties] AND alive[prop]'

Entrez.email = snakemake.params.email

# -- search for genes -- #
handle = Entrez.esearch(db="gene",
                        term=term,
                        retstart=0,
                        retmax=1)

# sleep as queries are limited to 3 a second
time.sleep(0.34)

record = Entrez.read(handle)

count = int(record["Count"])

ids = []
# NOTE: can not get more than 100.000 at once
for entry_ids in chunks(range(count), 100000):
    handle = Entrez.esearch(db="gene",
                            term=term,
                            retstart=entry_ids[0],
                            retmax=len(entry_ids))

    record = Entrez.read(handle)
    
    ids += record["IdList"]
    
    # sleep as queries are limited to 3 a second
    time.sleep(0.34)

# -- download genes -- #
# NOTE: split ids into chunks as fetching everything omits results
table_chunks = []
for id_chunk in chunks(ids, 1000):
    table_chunks.append(
        pd.read_csv(
            Entrez.efetch(db="gene", id=id_chunk, retmode="text", rettype="tabular"),
            sep="\t"
        )
    )
    
    # sleep as queries are limited to 3 a second
    time.sleep(0.34)

table = pd.concat(table_chunks)

# write table of protein coding genes
table.to_csv(snakemake.output.genes, sep="\t", index=False)
