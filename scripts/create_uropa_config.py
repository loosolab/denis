import json

# See https://uropa-manual.readthedocs.io/config.html
# For detailed explanation on config parameters

# load template
with open(snakemake.input.template, "r") as config_file:
    config = json.load(config_file)

# add input gtf and bed files
config["gtf"] = snakemake.input.gtf
config["bed"] = snakemake.input.bed

# write config as json
with open(snakemake.output[0], 'w') as outfile:
    json.dump(config, outfile, indent=4)
