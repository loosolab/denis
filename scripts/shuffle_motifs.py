# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import os
import random
from tobias.utils.motifs import MotifList, OneMotif

def make_MotifList(file_paths):
    '''
    Create a MotifList object from list of motif files.
    '''
    # load motifs
    motif_list = MotifList()
    
    for file in file_paths:
        motif_list = motif_list + MotifList().from_file(file)
    
    return MotifList(motif_list)

# create output dir if necessary
os.makedirs(snakemake.output[0], exist_ok=True)

# load motifs
motiflist = make_MotifList(snakemake.input)

# shuffle motif positions and bases
for motif in motiflist:
    counts = motif.counts
    # add pseudocount, otherwise positions could end up with all zero
    counts = [[j+1 for j in i] for i in counts]
    # double size by repeating motif once
    counts = [b * 2 for b in counts]
    # shuffle position
    [random.shuffle(b) for b in counts]
    # transpose matrix
    counts = list(map(list, zip(*counts)))
    # shuffle bases
    [random.shuffle(b) for b in counts]
    # transpose back
    counts = list(map(list, zip(*counts)))
    
    shuffled_motif = OneMotif(motifid=f"shuffled_{motif.id}", counts=counts, name=f"shuffled_{motif.name}")
    
    # save
    shuffled_motif.to_file(os.path.join(snakemake.output[0], f"{shuffled_motif.id}_{shuffled_motif.name}.meme"), fmt="meme")
    shuffled_motif.logo_to_file(os.path.join(snakemake.output[0], f"{shuffled_motif.id}_{shuffled_motif.name}.png"))
