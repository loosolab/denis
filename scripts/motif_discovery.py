# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import random
from kneed import KneeLocator
from scipy.optimize import curve_fit
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import subprocess
from timeit import default_timer as timer
from pathlib import Path
import re
from decimal import Decimal
from xml.etree import ElementTree
from lmfit import Model
import tqdm

from tobias.utils.motifs import MotifList

def run_meme(input_fasta, output_dir, threads=8, min_length=4, max_length=20, e_value=0.05, nmotifs=30, minsites=-1, maxsites=-1, mod="zoops", objfun="classic", log_file=None):
    '''
    Returns runtime in seconds.
    
    :param: input_fasta Input sites.
    :param: output_dir Meme output directory.
    :param: threads Number of parallel jobs.
    :param: min_length Minimum length of motifs.
    :param: max_length Maximum length of motifs.
    :param: e_value Stop when e_value above threshold.
    :param: nmotifs Stop when x motifs are found.
    :param: mod Distribution of motif sites. 'oops' = One Occurrence Per Sequence, 'zoops' = Zero or One Occurrence Per Sequence, 'anr' = Any Number of Repetitions. See https://meme-suite.org/meme/doc/meme.html?man_type=web
    :param: objfun How to calculate the motif's statistical significance, width and more. One of ['classic', 'de', 'se', 'cd', 'ce', 'nc']. see https://meme-suite.org/meme/doc/meme.html?man_type=web
    :param: minsites -1 for default or integer specifying minimum number of sites for motifs. see https://meme-suite.org/meme/doc/meme.html?man_type=web
    :param: maxsites -1 for default or integer specifying maximum number of sites for motifs. see https://meme-suite.org/meme/doc/meme.html?man_type=web
    :param: log_file If given will append all output to file.
    '''
    if not mod in ['oops', 'zoops', 'anr']:
        raise ValueError(f"Parameter mod has to be one of ['oops', 'zoops', 'anr']. Found {mod}")
    if not objfun in ['classic', 'de', 'se', 'cd', 'ce', 'nc']:
        raise ValueError(f"Parameter objfun has to be one of ['classic', 'de', 'se', 'cd', 'ce', 'nc']. Found {objfun}")
    
    meme_cmd = ["meme",
                input_fasta,
                "-o " + output_dir,
                "-dna",
                "-revcomp",
                f"-mod {mod}", # distribution of motif sites
                f"-objfun {objfun}",
                f"-minw {min_length}", # min motif length
                f"-maxw {max_length}", # max motif length
                f"-evt {e_value}", # stop when above e-value
                f"-nmotifs {nmotifs}", # stop when x motifs found
                # f"-wnsites 0", # bias towards number of specified sites
                f"-p {threads}", # number of threads
                f"-brief $(wc -l < {input_fasta})"] # include sequence list if below given value

    if maxsites > 0:
        meme_cmd.append(f"-maxsites {maxsites}")
    if minsites > 0:
        meme_cmd.append(f"-minsites {minsites}")

    # i'm too stupid to correctly use command list so join it
    meme_cmd = " ".join(meme_cmd)
    
    if log_file != None:
        # write command to log
        log_handle = open(log_file, 'a')
        
        log_handle.write("############################## MEME ##############################\n")
        log_handle.write("call:\n")
        log_handle.write(meme_cmd + "\n")
        log_handle.write("##################################################################\n")
        log_handle.flush() # To force write before subprocess output
    else:
        log_handle = None

    start = timer()
    subprocess.call(meme_cmd, shell=True, stdout=log_handle, stderr=log_handle)
    end = timer()
    
    if log_file != None:
        log_handle.close()
    
    return end - start

def make_MotifList(file_paths):
    '''
    Create a MotifList object from list of motif files.
    '''
    # load motifs
    motif_list = MotifList()
    
    for file in file_paths:
        motif_list = motif_list + MotifList().from_file(file)
    
    return MotifList(motif_list)

def distance_matrix(file_paths, output_file):
    '''
    Create and save a all vs. all distance matrix.
    '''
    # create motif_list
    motifs = make_MotifList(file_paths)
        
    # cluster then save similarity matrix
    motifs.cluster(threshold=0.3, metric="pcc", clust_method="average")
        
    # create path if needed
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    motifs.similarity_matrix.to_csv(output_file, sep = '\t')

def distance_heatmap_all(distance_matrix, output_file="all_heatmap.png"):
    '''
    A heatmap with all vs. all.
    '''
    
    # load distance table
    dist_matrix = pd.read_csv(distance_matrix, sep="\t", index_col=0)
        
    # plot full heatmap
    plot = sns.clustermap(dist_matrix,
                        cmap="YlOrRd_r",
                        vmin=0,
                        vmax=1#,
                        #figsize=(x,y)
                        )
    plot.ax_heatmap.set_xlabel("All Motifs")
    plot.ax_heatmap.set_ylabel("All Motifs")
    plot.savefig(output_file, bbox_inches='tight', dpi=100)
    plt.close()

def extract_motif_info(file):
    '''
    Extract all information from the headlines of each motif in the given meme.txt.
    
    Output is a dict where the motif names are keys, containing a dict with the information.
    Keys of the second dict are: width, sites, llr (loglikelihood ratio), E-value
    '''
    info = dict()

    with open(file, "r") as f:
        for line in f:
            if "llr" in line:
                motif_name = re.search("^MOTIF (.*)\t", line).group(1)

                # remove leading name
                single_motif_info = line.split("\t")[1]
                # split into entries
                single_motif_info = re.split("(\S+ =\s+\S+)\s+", single_motif_info)
                # remove empty strings
                single_motif_info = list(filter(None, single_motif_info))
                # create dict
                single_motif_info = dict(re.split("\s+=\s+", x) for x in single_motif_info)

                info[motif_name] = single_motif_info
    
    return info

def create_motif_stats(meme_file):
    '''
    Create and return a matrix of values for each motif.
    nsites, nsites percentage, distance, e-value
    '''
    motifs = MotifList().from_file(meme_file)
    
    motif_stats = pd.DataFrame(columns=["id", "name"])
    
    addon_info = extract_motif_info(meme_file)
    
    # create info matrix row by row
    for motif in motifs:
        motif.information_content()
        motif.gc_content()
        
        motif_info = pd.DataFrame(
            {'id': [motif.id],
             'name': motif.name, 
             'information_content': motif.ic,
             'GC%': motif.gc,
             'nsites': int(motif.info['nsites']),
             'e_value': float(motif.info['E']),
             'log_likelihood': int(addon_info[f"{motif.id} {motif.name}"]['llr']),
             'width': int(motif.info['w'])}
        )

        motif_stats = pd.concat([motif_stats, motif_info], ignore_index=True)
    
    # add percentage
    motif_stats["nsites percentage"] = motif_stats["nsites"] / motif_stats["nsites"].sum()
    
    return motif_stats

def knee_plot(motif_stats, on, output_file="knee_plot.png", log=True, curve="concave", direction="increasing", fit_function="linear"):
    '''
    Create a knee plot and return the knee-value threshold to filter on.
    #TODO reimplement parameters log, curve and direction
    '''
    motif_stats = motif_stats.copy()
    
    motif_stats.sort_values(by=on, ascending=True, inplace=True)
    motif_stats[on + "_order"] = range(1, len(motif_stats)+1)

    if log:
        motif_stats["log_" + on] = np.log10(motif_stats[on])
        
        # replace '-inf' with -400 as python seems to lack the precision
        motif_stats.loc[motif_stats["log_" + on] <= -np.inf, "log_" + on] = -400
    
    # fit curve
    # see https://stackoverflow.com/a/42793748    
    # define fit function and create model
    if fit_function == "linear":
        def fit_func(x, m1, b1, m2, b2):
            # get x where lines intersect
            m = m2 - m1
            b = b1 - b2
            intersect_x = np.nan if m == 0 else b / m
            
            result = []
            
            for xi in x:
                if intersect_x is np.nan or intersect_x >= xi:
                    result.append(m1 * xi + b1)
                else:
                    result.append(m2 * xi + b2)
            
            return result
        
        model = Model(fit_func)
        params = model.make_params(m1=-np.min(motif_stats["log_" + on])/(len(motif_stats["log_" + on])/2),
                                   b1=np.min(motif_stats["log_" + on]),
                                   m2=0,
                                   b2=0)
        
    elif fit_function == "curved":
        if log:
            def fit_func(x, a, b, c):
                return a * np.log(b * x) + c
        else:
            def fit_func(x, a, b, c):
                return np.exp((x - c ) / a) / b
        
        # TODO set default params for not log
        model = Model(fit_func)
        params = model.make_params(a=1, c=-1)
        params.add(name='b', min=0.1)
            
    # fit model
    result = model.fit(motif_stats["log_" + on], params=params, x=motif_stats[on + "_order"])

    #if log:
    #    popt, pcov = curve_fit(fit_func, motif_stats[on + "_order"], motif_stats["log_" + on], bounds=([-np.inf, 0, -np.inf], np.inf))
    #else:
    #    popt, pcov = curve_fit(fit_func, motif_stats[on + "_order"], motif_stats[on], bounds=([-np.inf, 0, -np.inf], np.inf))
    
    # see https://stackoverflow.com/questions/51762514/find-the-elbow-point-on-an-optimization-curve-with-python
    kn_smooth = KneeLocator(motif_stats[on + "_order"], result.best_fit, curve=curve, direction=direction)
    smooth_threshold = kn_smooth.knee

    # create knee plot
    plt.figure()
    plt.ylabel(on)
    
    plt.plot(motif_stats[on + "_order"], result.init_fit, 'r--', label='Initial fit')
    plt.plot(motif_stats[on + "_order"], result.best_fit, 'r-', label="Best fit")
    if smooth_threshold != None:
        plt.vlines(smooth_threshold, plt.ylim()[0], plt.ylim()[1], linestyles='dotted', colors='tab:orange', label='Knee')
    
    if log:
        plt.plot(motif_stats[on + "_order"], motif_stats["log_" + on], 'bo')
    else:
        plt.plot(motif_stats[on + "_order"], motif_stats[on], 'bo')
    # plt.vlines(threshold, plt.ylim()[0], plt.ylim()[1], linestyles='dashed', colors='b')
    plt.legend(loc='best')
    plt.savefig(output_file, dpi=800)
    plt.close()
    
    if smooth_threshold == None:
        return None
    
    return motif_stats.loc[motif_stats[on + "_order"] == smooth_threshold, on].values[0]

def invert_dict(d):
    '''
    Switch elements with keys and return.
    '''
    reverse = {}
    for key in d:
        for value in d[key]:
            reverse[value] = reverse.get(value, []) + [key]
        
    return reverse

def filter_sites(xml_file,
                 fasta_file,
                 motif_names,
                 output_file="filtered.fasta",
                 motif_output_dir=None,
                 save="all",
                 prefix="",
                 accepted=True,
                 random_p=0.1,
                 reuse_sites=False,
                 min_len=4):
    '''
    Remove sites from fasta that are in the xml.
    
    :param: xml_file Path to xml file containing motif to sequence relation
    :param: fasta_file Path to input fasta. All sequences available to meme.
    :param: motif_names Name of motif which sites should be removed. If empty list removes given percentage (See param random_p) of sites at random.
    :param: output_file Fasta without the sequences that where used in a motif (occured in xml)
    :param: motif_output_dir If set write fasta for each motif.
    :param: save If "all" all will be saved. If "accepted" only accepted motifs (in motif_names) will be stored. If 'none' no motifs will be stored.
    :param: prefix Prefix to files stored in motif_output_dir
    :param: accepted If True add prefix "accepted" to all files of motifs declared in motif_names
    :param: random_p Percentage of sites that will be randomly removed if motif_names length is 0. Default = 0.1 (10%).
    :param: reuse_sites Boolean, if True will split sites at motif match location to reuse parts. If False sequence will be discarded.
    :param: min_len int; Minimum length a sequence can have to be stored if reuse_sites=True.
    '''
    # get number of records in fasta
    fasta_len = len([1 for line in open(fasta_file) if line.startswith(">")])

    # remove sites randomly
    if len(motif_names) < 1:
        # index of sites to remove
        remove_index = random.sample(range(fasta_len), fasta_len * random_p)
        
        # load fasta
        file = SeqIO.parse(fasta_file, "fasta")
        
        # filter and write fasta
        with open(output_file, "w") as out:
            for index, record in tqdm.tqdm(enumerate(file), total=fasta_len):
                if index in remove_index:
                    continue
                
                SeqIO.write(record, out, "fasta")

    # remove used sequences
    else:
        ##### load used sequences #####
        # load xml
        xml = ElementTree.parse(xml_file)
        root = xml.getroot()

        # fasta sequence id -> xml sequence id
        fasta_xml = {seq.attrib["name"]: seq.attrib["id"] for seq in root.find("training_set").findall("sequence")}
        # xml sequence id -> fasta sequence id
        xml_fasta = {seq.attrib["id"]: seq.attrib["name"] for seq in root.find("training_set").findall("sequence")}
        # fasta sequence id -> [motif binding sequence(s)]
        motif_bs = {}
        # motif name -> [xml sequence id]; sequences used by a motif
        motif_xml = {}
        # fasta seq ids of records belonging to an accepted motif; will be removed
        omit_rec = set()
        

        # populate motif_bs, motif
        for motif in root.find("motifs").findall("motif"):
            name = motif.attrib["name"] + " " + motif.attrib["alt"]

            for site in motif.find("contributing_sites").findall("contributing_site"):
                site_id = site.attrib["sequence_id"]

                # get motif associated sequence ids
                motif_xml.setdefault(name, []).append(site_id)

                # get motif binding sequence
                seq = "".join([letter.attrib["letter_id"] for letter in site.find("site")])

                if site.attrib["strand"] == "minus":
                    seq = str(Seq(seq).reverse_complement())

                motif_bs.setdefault(xml_fasta[site_id], []).append(seq)

                if name in motif_names:
                    omit_rec.add(xml_fasta[site_id])

        ##### fasta #####
        fasta = SeqIO.parse(fasta_file, "fasta")

        # individual motif fastas
        # write all files add once for better performance
        ind_fastas = {}
        full_fasta = []

        # iterate over all records in fasta
        for record in tqdm.tqdm(fasta, total=fasta_len):
            # write individual motif fasta
            if motif_output_dir:
                # create motif output dir if needed
                os.makedirs(motif_output_dir, exist_ok=True)

                for name, seq_ids in motif_xml.items():
                    # Skip if seq does not belong to an accepted motif and save all is off or none should be saved.
                    if save in ["accepted"] and not name in motif_names or save in ["none"]:
                        continue

                    # skip sequences not associated with the motif
                    if fasta_xml[record.id] not in seq_ids:
                        continue

                    # add "accepted" prefix to individual motif output file
                    tmp_prefix = prefix
                    if name in motif_names and accepted:
                        tmp_prefix = prefix + "accepted_"

                    # create motif fasta + path and remove whitespaces
                    motif_fasta = os.path.join(motif_output_dir, tmp_prefix + name + ".fasta")
                    motif_fasta = motif_fasta.replace(" ", "_")

                    # store seq for later writing
                    ind_fastas.setdefault(motif_fasta, []).append(record)

            # write to full output fasta
            if reuse_sites and record.id in omit_rec:
                # remove motif binding location and keep flanking sequences
                for bind_seq in motif_bs[record.id]:
                    # split sequence
                    left_flank, right_flank = record.seq.split(bind_seq, 1)

                    # prepare ids; assumes fasta id format "[chrX]:[start]-[end]"
                    id_list = record.id.split(":")
                    start, end = map(int, id_list[-1].split("-"))
                    id_left = ":".join(id_list[:-1] + [f"{start}-{start + len(left_flank)}"])
                    id_right = ":".join(id_list[:-1] + [f"{end - len(right_flank)}-{end}"])

                    # generate record for each flank
                    if len(left_flank) >= min_len:
                        l_rec = SeqRecord(
                            seq=Seq(left_flank),
                            id=id_left,
                            name=id_left,
                            description=id_left
                        )

                        full_fasta.append(l_rec)

                    if len(right_flank) >= min_len:
                        r_rec = SeqRecord(
                            seq=Seq(right_flank),
                            id=id_right,
                            name=id_right,
                            description=id_right
                        )

                        full_fasta.append(r_rec)
            elif record.id not in omit_rec:
                # write full record
                full_fasta.append(record)

        # write full fasta
        with open(output_file, "w") as out:
            SeqIO.write(full_fasta, out, "fasta")

        # write individual motif fastas
        for m_fasta, records in ind_fastas.items():
            with open(m_fasta, "w") as out_file:
                SeqIO.write(records, out_file, "fasta")

def saveMotifs(input_file, output_dir, motif_names, prefix):
    '''
    Save motifs from given motif files into seperate meme files.
    
    :param: input_file .meme or .txt of motifs
    :param: output_dir Where the files should be saved.
    :param: motif_names Which motifs of input should be saved. List or string
    :param: prefix to add before "[motif_name].meme"
    '''
    # cast to list if needed
    motif_names = motif_names if isinstance(motif_names, list) else [motif_names]
    
    # load file
    motif_list = MotifList().from_file(input_file)
    
    for motif in motif_list:
        name = f"{motif.id} {motif.name}"
        
        if name in motif_names:
            output_file = os.path.join(output_dir, f"{prefix}_{name}.meme")
            output_file = output_file.replace(" ", "_")
            
            motif.to_file(output_file, fmt="meme")

def motif_discovery(fasta, 
                    output_dir, 
                    retry=3, 
                    escape_method='strongest', 
                    random_sites=0.1, 
                    save_discarded=False, 
                    threads=1, 
                    meme_log=None, 
                    filter_var="e_value", 
                    fit_func="linear", 
                    min_len=4, 
                    max_len=20, 
                    threshold=0.05, 
                    mod="zoops", 
                    objfun="classic",
                    nmotifs=30,
                    minsites=-1,
                    maxsites=-1,
                    maxiter=-1,
                    reuse_sites=True):
    '''
    Run motif discovery
    
    Iterative process of motif discovery.
    Run Meme and filter with a knee plot.
    Sites of accepted motifs will be removed prior to next iteration.
    
    Will stop when:
        - no motifs are found
        - no sites are left
        - no filter threshold (knee) for 'retry' iterations
        - no motifs are accepted
    
    :param fasta: Input of fasta file used in first iteration.
    :param output_dir: Where all results will be saved.
    :param retry: Number of iterations with no filter (knee) until discovery is stopped.
    :param escape_method: one of ['strongest', 'weakest', 'random']. How to escape a local optimum. Either remove the sites corresponding to the strongest/ weakest motif or remove 10% random sites.
    :param random_sites: If escape_method is 'random' remove this percentage of sites. Default = 0.1 (10%).
    :param save_discarded: Whether the not accepted motifs should be saved along side the accepted ones.
    :param threads: Number of threads to use when running meme.
    :param meme_log: Redirect meme output to given file.
    :param filter_var: One which variable should be filtered. Either "e_value" or "log_likelihood". #TODO implement log_likelihood
    :param fit_func: Either ["linear", "curved"]. Decides the shape of function used for fitting and subsequent knee finding.
    :param min_len: Minimum width of motifs.
    :param max_len: Maximum width of motifs.
    :param threshold: Threshold at which MEME will stop searching for motifs.
    :param mod: Distribution of motif sites. 'oops' = One Occurrence Per Sequence, 'zoops' = Zero or One Occurrence Per Sequence, 'anr' = Any Number of Repetitions. See https://meme-suite.org/meme/doc/meme.html?man_type=web
    :param objfun: How to calculate the motif's statistical significance, width and more. One of ['classic', 'de', 'se', 'cd', 'ce', 'nc']. see https://meme-suite.org/meme/doc/meme.html?man_type=web
    :param nmotifs: Stop when x motifs are found.
    :param minsites: -1 for default or integer specifying minimum number of sites for motifs. see https://meme-suite.org/meme/doc/meme.html?man_type=web
    :param maxsites: -1 for default or integer specifying maximum number of sites for motifs. see https://meme-suite.org/meme/doc/meme.html?man_type=web
    :param maxiter: Stop after x iterations. -1 = unlimited
    :param: reuse_sites Boolean, if True will split sites at motif match location to reuse parts. If False sequence will be discarded.
    '''

    finished = False
    attempts = 0

    # number of current iteration
    iteration = 0

    start_time = timer()

    while not finished:
        if maxiter > 0 and iteration > maxiter:
            print("Stopped because maximum iteration reached.")
            finished = True
            break
        
        print("Iteration: " + str(iteration))
        output_it = os.path.join(output_dir, f"iteration_{iteration}")
        os.makedirs(output_it, exist_ok=True)
        
        ##### RUN MEME #####
        print("Meme:")
        meme_time = run_meme(fasta, os.path.join(output_it, "meme"), threads=threads, min_length=min_len, max_length=max_len, e_value=threshold, mod=mod, objfun=objfun, nmotifs=nmotifs, minsites=minsites, maxsites=maxsites, log_file=meme_log)
        print(str(meme_time / 60) + " min.")
        
        # stop when no motifs are found
        if not glob.glob(os.path.join(output_it, "meme", "*.png")):
            print("Stopped because no motifs were found.")
            finished = True
            break

        ##### MOTIF STATISTICS #####
        print("Make motif statistics")
        os.makedirs(os.path.join(output_it, "validation"), exist_ok=True)
        
        stat_matrix = create_motif_stats(os.path.join(output_it, "meme", "meme.txt"))
        
        ##### FILTER MOTIFS #####
        print("Filter motifs:")
        
        # keep going when there is only one result
        # stat_matrix has to have at least one more datapoint than fitting function has independent vars
        if fit_func == "linear" and len(stat_matrix) > 3 or fit_func == "curved" and len(stat_matrix) > 2:
            knee_threshold = knee_plot(stat_matrix,
                                    on=filter_var,
                                    output_file=os.path.join(output_it, "validation", "knee_plot.png"),
                                    log=False if filter_var == "log_likelihood" else True,
                                    curve="convex" if filter_var == "log_likelihood" else "concave",
                                    direction="increasing",
                                    fit_function=fit_func)
        else:
            knee_threshold = None

        # if no knee found stop after 3 iteration
        if knee_threshold is None:
            if attempts >= retry:
                finished = True
                knee_threshold = 0 if filter_var == "log_likelihood" else 1
            else:
                attempts += 1
                if escape_method == "strongest":
                    knee_threshold = stat_matrix["log_likelihood"].nlargest(1).iloc[0] if filter_var == "log_likelihood" else stat_matrix["e_value"].nsmallest(1).iloc[0]
                elif escape_method == "weakest":
                    knee_threshold = stat_matrix["log_likelihood"].nsmallest(1).iloc[0] if filter_var == "log_likelihood" else stat_matrix["e_value"].nlargest(1).iloc[0]
        else:
            attempts = 0
        
        # add accepted motifs column
        stat_matrix['accepted'] = (stat_matrix["log_likelihood"] >= knee_threshold) if filter_var == "log_likelihood" else (stat_matrix["e_value"] <= knee_threshold)
        
        # get accepted motifs
        accepted_motifs = stat_matrix[stat_matrix['accepted']]
        accepted_names = (accepted_motifs["id"] + " " + accepted_motifs["name"]).tolist()
        
        # save matrix
        stat_matrix.to_csv(os.path.join(output_it, "validation", "motif_stats.tsv"), sep="\t", index=False)
        
        # if there is nothing accepted stop
        if len(accepted_motifs) < 1 and not escape_method == "random":
            finished = True
        
        # filter sites for next iteration and save used sites
        filter_time = timer()
        filter_sites(xml_file=os.path.join(output_it, "meme", "meme.xml"),
                    fasta_file=fasta,
                    motif_names=accepted_names,
                    output_file=os.path.join(output_it, "meme", "filtered.fasta"),
                    motif_output_dir=os.path.join(output_dir, "1_motifs"),
                    save="all" if save_discarded else "accepted" if attempts == 0 else "none",
                    prefix="iteration_" + str(iteration) + "_",
                    accepted=attempts == 0,
                    random_p=random_sites,
                    reuse_sites=reuse_sites,
                    min_len=min_len)
        
        # set accepted False if attempts > 0
        if attempts > 0:
            stat_matrix['accepted'] = False
            accepted_names = []

        print(str((timer() - filter_time) / 60) + " min.")
        
        ##### SAVE MOTIFS #####
        print("Write motif files")
        
        # save accepted
        if attempts == 0:
            saveMotifs(input_file=os.path.join(output_it, "meme", "meme.txt"),
                    output_dir=os.path.join(output_dir, "1_motifs"),
                    motif_names=accepted_names, 
                    prefix="iteration_" + str(iteration) + "_accepted")

        # save not accepted
        if save_discarded:
            not_accepted = list(set((stat_matrix["id"] + " " + stat_matrix["name"]).tolist()) - set(accepted_names))
            saveMotifs(input_file=os.path.join(output_it, "meme", "meme.txt"),
                    output_dir=os.path.join(output_dir, "1_motifs"),
                    motif_names=not_accepted, 
                    prefix="iteration_" + str(iteration))
        
        # setup next run
        fasta = os.path.join(output_it, "meme" , "filtered.fasta")
        iteration += 1

    if not os.path.isdir(os.path.join(output_dir, "1_motifs")) or len(next(os.walk(os.path.join(output_dir, "1_motifs")))[2]) <= 0:
        raise Exception("No motifs discovered!")

    print("Total time: " + str((timer() - start_time) / 60) + " min.")

if __name__ == "__main__":
    # parse command line args
    import argparse
    import sys
    
    parser = argparse.ArgumentParser(description='Find motifs out of fasta sequences. Iterative process of finding motifs and removing corresponding sites.')
    
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('-f', '--fasta', required=True, help='The input fasta file.')
    required.add_argument('-o', '--output_dir', required=True, help='Path to the output directory.')
    optional.add_argument('-r', '--retry', default=3, type=int, help='Number of times to retry when no knee value is found.')
    optional.add_argument('-e', '--escape_method', default='strongest', choices=['strongest', 'weakest', 'random'], help='Choose which method to use when no knee value as a filter threshold is found.')
    optional.add_argument('-s', '--save_discarded', action='store_true', help='If present save discarded motifs.')
    optional.add_argument('-t', '--threads', default=1, type=int, help='Number of threads to use.')
    optional.add_argument('-l', '--meme_log', default=None, help='Will redirect meme output to given file.')
    optional.add_argument('-k', '--fit_func', default='curved', choices=['linear', 'curved'], help='What equation should be used for fitting and knee finding.')
    optional.add_argument('-n', '--min_len', default=4, type=int, help='Minimum motif width.')
    optional.add_argument('-m', '--max_len', default=20, type=int, help='Maximum motif width.')
    optional.add_argument('-v', '--threshold', default=0.05, type=float, help='E-value threshold at which MEME will stop searching for motifs.')
    optional.add_argument('-d', '--mod', default='zoops', choices=['oops', 'zoops', 'anr'], help="Distribution of motif sites. 'oops' = One Occurrence Per Sequence, 'zoops' = Zero or One Occurrence Per Sequence, 'anr' = Any Number of Repetitions. See https://meme-suite.org/meme/doc/meme.html?man_type=web")
    optional.add_argument('-b', '--objfun', default='classic', choices=['classic', 'de', 'se', 'cd', 'ce', 'nc'], help="How to calculate the motif's statistical significance, width and more. see https://meme-suite.org/meme/doc/meme.html?man_type=web")
    optional.add_argument('-a', '--nmotifs', default=30, type=int, help='Stop when x motifs are found.')
    optional.add_argument('-c', '--minsites', default=-1, type=int, help='Integer specifying minimum number of sites for motifs. see https://meme-suite.org/meme/doc/meme.html?man_type=web. -1 for default.')
    optional.add_argument('-g', '--maxsites', default=-1, type=int, help='Integer specifying maximum number of sites for motifs. see https://meme-suite.org/meme/doc/meme.html?man_type=web. -1 for default.')
    optional.add_argument('-i', '--maxiter', default=-1, type=int, help='Stop after x iterations. -1 = unlimited')

    args = parser.parse_args()
    
    # print help if no parameters given
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    
    args = vars(args)

    # run
    motif_discovery(**args)
