# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import random
import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tobias.utils.motifs import OneMotif

def random_motif(fasta, output_dir, prefix="", quantity=5, min_len=4, max_len=20, percentage=0.2, multiple=False):
    """
    Create a number of random motifs from fasta.
    
    Parameters:
        fasta (string): Path to fasta file.
        output_dir (string): Path to output directory.
        prefix (string): Prefix that will be added to all file names.
        quantity (int): Number of motifs to generate.
        min_len (int): Minimum motif size.
        max_len (int): Maximum motif size.
        percentage (float): Percentage of fasta entries used per motif. For multiple=True percentage of bases used per motif.
        multiple (boolean): If 'True' will allow multiple sites per fasta entry.
    """
    # create output dir if necessary
    os.makedirs(output_dir, exist_ok=True)
    
    step = (max_len - min_len) / (quantity - 1)
    sizes = [round(min_len + i * step) for i in range(quantity)]
    
    for i, size in enumerate(sizes):
        id = f"{prefix}_{i}"
        name = f"random"
        f_basename = f"{id}_{name}"

        # write motif fasta and bed
        out_fasta = os.path.join(output_dir, f"{f_basename}.fasta")
        with open(out_fasta, "w") as out:
            with open(os.path.join(output_dir, f"{f_basename}.bed"), "w") as out_bed:
                for record in SeqIO.parse(fasta, "fasta"):
                    if not multiple:
                        # single site per fasta entry
                        if "N" in record.seq or size > len(record.seq) or random.random() > percentage:
                            continue
                        
                        # select random site
                        x = random.randint(0, len(record.seq) - size)
                        record.seq = record.seq[x:x+size]
                        
                        # adjust id
                        match = re.search(r"(\d+)-(\d+)$", record.id)
                        record.id = re.sub(r"(\d+)-(\d+)$", f"{int(match.group(1)) + x}-{int(match.group(1)) + x+size}", record.id)
                        record.name = record.id
                        record.description = record.id
                        
                        # write fasta
                        SeqIO.write(record, out, "fasta")
                        
                        # write bed
                        bed_match = re.search(r"(\w+):(\d+)-(\d+)$", record.id)
                        out_bed.write("\t".join(bed_match.groups()) + "\n")
                    else:
                        # multiple sites per fasta entry
                        num_sites = round(len(record.seq) * percentage / size)
                        
                        # radomly select sites until enough
                        saved = 0
                        while saved < num_sites:
                            # select random site
                            x = random.randint(0, len(record.seq) - size)
                            
                            if "N" in record.seq[x:x+size]:
                                continue
                            
                            # create fasta subset entry
                            info = f"{record.id}:{x}-{x+size}"
                            subset = SeqRecord(seq=record.seq[x:x+size], id=info, name=info, description=info)
                            
                            # write fasta
                            SeqIO.write(subset, out, "fasta")
                            
                            # write bed
                            out_bed.write("\t".join(map(str, [record.id, x, x+size])) + "\n")
                            
                            saved += 1
            
        # write motif meme
        motif = OneMotif.from_fasta(out_fasta, motifid=id, name=name)
        
        motif.logo_to_file(os.path.join(output_dir, f"{f_basename}.png"), ylim="auto")
        
        motif.to_file(os.path.join(output_dir, f"{f_basename}.meme"), fmt="meme")


if __name__ == "__main__":
    # parse command line args
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(description='Generate an arbitrary number of random motifs from a fasta.')
    
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('-f', '--fasta', required=True, help='The input fasta file.')
    required.add_argument('-o', '--output_dir', required=True, help='Path to the output directory.')
    optional.add_argument('-r', '--prefix', default='', help='Prefix that will be added to all file names.')
    optional.add_argument('-q', '--quantity', default=5, type=int, help='Number of motifs to generate.')
    optional.add_argument('-i', '--min_len', default=4, type=int, help='Minimum motif size.')
    optional.add_argument('-a', '--max_len', default=20, type=int, help='Maximum motif size.')
    optional.add_argument('-p', '--percentage', default=0.2, type=float, help='Percentage of fasta entries used per motif. For multiple=True percentage of bases used per motif')
    optional.add_argument('-m', '--multiple', action='store_true', help='If present allow for multiple sites per fasta entry.')

    args = parser.parse_args()
    
    # print help if no parameters given
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    
    args = vars(args)

    # run
    random_motif(**args)