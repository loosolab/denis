import os
import glob
import pandas as pd

from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.godag_plot import plot_results


def ensmbl_to_entrez(ensmbl_list, mapping_table):
    """ Convert list of Ensembl GeneID to Enrez GeneID. Multi-matches are kept and NA ids omitted. """  
    entrez_list = []
    for ens_id in ensmbl_list:
        val = list(mapping_table.loc[mapping_table["Ensembl_gene_identifier"] == ens_id, "GeneID"].values)

        if val:
            entrez_list.extend(val)
        
    return entrez_list


def goea_table(go_list):
    """ Format list of go-terms into a pandas DataFrame. """
    ns = {"BP": "biological_process", "MF": "molecular_function", "CC": "cellular_component"}
    en = {"e": "enriched", "p": "purified"}

    table = {}

    for go_term in go_list:
        fold_enr = (go_term.study_count / go_term.study_n) / (go_term.pop_count / go_term.pop_n)

        table.setdefault("GO", []).append(go_term.GO)
        table.setdefault("name", []).append(go_term.name)
        table.setdefault("namespace", []).append(ns[go_term.NS])
        table.setdefault("depth", []).append(go_term.depth)
        table.setdefault("enrichment", []).append(en[go_term.enrichment])
        table.setdefault("study_count", []).append(go_term.study_count)
        table.setdefault("study_n", []).append(go_term.study_n)
        table.setdefault("pop_count", []).append(go_term.pop_count)
        table.setdefault("pop_n", []).append(go_term.pop_n)
        table.setdefault("p-value", []).append(go_term.p_uncorrected)
        table.setdefault("fold enrichment", []).append(fold_enr)
        table.setdefault("FDR", []).append(go_term.p_fdr_bh)
        table.setdefault("GO_alt_ids", []).append(go_term.goterm.alt_ids)

    return pd.DataFrame(table)

# ---------- Setup ---------- #

# load UROPA annotation
annotations = pd.read_csv(glob.glob(os.path.join(snakemake.input.annotation, "*_allhits.txt"))[0], sep="\t")

# load background genes
bg_genes = pd.read_csv(snakemake.input.bg_genes, sep="\t")

# load gene2ensembl for entrez to ensembl id mapping
gene2ensembl = pd.read_csv(snakemake.input.gene2ensembl, sep="\t")
# subset to needed organism (for speed up)
gene2ensembl = gene2ensembl[gene2ensembl["#tax_id"] == bg_genes["tax_id"].values[0]]

# parse & load go-term associations
obo_graph = GODag(snakemake.input.obo, prt=None)

# parse & load gene to go-term annotation
gene2go_annot = Gene2GoReader(snakemake.input.gene2go, taxids=[bg_genes["tax_id"].values[0]], prt=None)

# Get namespace2association where:
#    namespace is:
#        BP: biological_process
#        MF: molecular_function
#        CC: cellular_component
#    assocation is a dict:
#        key: NCBI GeneID
#        value: A set of GO IDs associated with that gene
ns2assoc = gene2go_annot.get_ns2assc()

# ---------- Gene Ontology Enrichment Analysis ---------- #

# setup enrichment analysis
goeaobj = GOEnrichmentStudyNS(
        bg_genes["GeneID"], # List of protein-coding genes
        ns2assoc, # geneid/GO associations
        obo_graph, # Ontologies
        propagate_counts = False,
        alpha=0.05, # default significance cut-off
        methods=['fdr_bh'], # use false discovery rate benjamini-hochberg procedure
        prt=None)

# run enrichment analysis & find enriched GO terms
go_terms = goeaobj.run_study(
    ensmbl_to_entrez(annotations["gene_id"].tolist(), gene2ensembl),
    prt=None # make silent
)

# filter significant
significant_go_terms = [r for r in go_terms if r.p_fdr_bh < 0.05]

# create dir if needed
os.makedirs(snakemake.output.dir, exist_ok=True) 

# save significant table
table = goea_table(significant_go_terms)
table.to_csv(os.path.join(snakemake.output.dir, "significant_go_terms.tsv"), sep="\t", index=None)

# save go network
plot_results(fout_png=os.path.join(snakemake.output.dir, "graph_{NS}.png"), goea_results=significant_go_terms, prt=None)
