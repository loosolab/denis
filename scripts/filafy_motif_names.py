from tobias.utils.utilities import filafy
from tobias.utils.motifs import MotifList
import pandas as pd

# load motifs
motifs = MotifList().from_file(snakemake.input[0])

table = {}

for m in motifs:
    table.setdefault("id", []).append(m.id)
    table.setdefault("name", []).append(m.name)
    table.setdefault("f_id", []).append(filafy(m.id))
    table.setdefault("f_name", []).append(filafy(m.name))

table = pd.DataFrame(table)

table["id_name"] = table["id"] + "_" + table["name"]
table["filafy_id_name"] = table["f_id"] + "_" + table["f_name"]

table.to_csv(snakemake.output[0], sep="\t", index=False)
