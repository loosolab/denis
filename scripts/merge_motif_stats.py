import pandas as pd
import re

table_list = []
for file in snakemake.input.stat_files:
    table = pd.read_csv(file, sep="\t")
    
    table['iteration'] = re.search(r'iteration_(\d+)', file).group(1)
    
    table_list.append(table)
    
merged = pd.concat(table_list)

merged.to_csv(snakemake.output[0], sep="\t", index=False)
