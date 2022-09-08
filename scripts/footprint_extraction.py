
import os.path
import time
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline

import pyBigWig
from Bio import SeqIO
from tobias.utils.motifs import MotifList, OneRegion

# fix 'Could not connect to any X display.'
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
import seaborn as sns


def plot_footprints(scores, slope, percentile_threshold, footprints, title, motifs=None):
    """
    Creates figure with two subplots. 
    Left plot = peak region with found footprints and motifs.
    Right plot = slope distribution for peak extensions.

    Parameters
    ----------
    scores : list
        Array of score per position.
    slope : list
        Array of slope per position.
    percentile_threshold : float
        Accepted slope range beginning from 50% Quantile. E.g. for 35 all slopes between the 15% and 85% quantile will be accepted.
    footprints : dict
        Dictionary with keys 'start' and 'end' which are lists of footprint positions.
    title : str
        Plot title.
    motifs : dict
        Dictionary same structure as in footprints parameter.

    Returns
    -------
    matplotlib.figure.Figure
        Fig object
    """
    fig, axs = plt.subplots(ncols=2, figsize=(14, 5), tight_layout={'rect': (0, 0, 1, 0.95)}) # prevent label clipping; leave space for title
    fig.suptitle(title, fontsize=16)

    # footprint plot
    axs[0].set_title('Extracted footprints (green)')
    color = 'tab:blue'
    axs[0].patch.set_alpha(0)
    axs[0].zorder = 10
    axs[0].set_xlabel('basepairs')
    axs[0].set_ylabel('binding score', color=color)
    axs[0].plot(scores, color=color)
    axs[0].tick_params(axis='y', labelcolor=color)
    # add peak rectangles
    for left, right in zip(footprints["start"], footprints["end"]):#zip(properties["left_edges_extended"], properties["right_edges_extended"]):
        # set linewidth small or else rectangles will overlap
        axs[0].add_patch(Rectangle((left, 0), width=right + 1 - left, height=np.max(scores)* 1.1, alpha=0.2, facecolor="tab:green", edgecolor='black', linewidth=1))
    # add motif rectangles
    if motifs is not None:
        for left, right in zip(motifs["start"], motifs["end"]):#zip(properties["motif_left"], properties["motif_right"]):
            # set linewidth small or else rectangles will overlap
            axs[0].add_patch(Rectangle((left, 0), width=right + 1 - left, height=np.max(scores)* 1.1, alpha=0.2, facecolor="tab:red", edgecolor='black', linewidth=1))
        
    # add second y-axis
    ax2 = axs[0].twinx()

    color2 = 'tab:orange'
    ax2.zorder = 5
    ax2.set_ylabel('1st order derivative', color=color2)
    ax2.plot(slope, color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)

    # distribution plot
    axs[1].set_title('First order derivative distribution')
    sns.histplot(slope, ax=axs[1], kde=True)
    # threshold lines
    lt = np.percentile(slope, 50 - percentile_threshold)
    axs[1].axvline(x=lt, label=f"{50 - percentile_threshold}% percentile = {np.round(lt, 10)}", color='red', linestyle='--')
    axs[1].axvline(x=np.median(slope), label=f"Median = {np.round(np.median(slope), 10)}", color='grey', linestyle='--')
    ut = np.percentile(slope, 50 + percentile_threshold)
    axs[1].axvline(x=ut, label=f"{50 + percentile_threshold}% percentile = {np.round(ut, 10)}", color='orange', linestyle='--')

    axs[1].legend() #show legend

    return fig


def extend_peak(scores, slopes, peaks, min_slope_threshold, max_slope_threshold, max_gap, gap_depth, reverse=False):
    """
    Extend footprint peaks on one side.

    Parameters
    ----------
    scores : list
        Array of score per position.
    slopes : list
        Array of slope per position.
    peaks : list
        Array of peak edge positions.
    min_slope_threshold : float
        Minimum accepted slope for extension.
    max_slope_threshold : float
        Maximum accepted slope for extension.
    max_gap maximum : int
        Accepted gap width to merge neighboring peaks.
    gap_depth : float
        Maximum accepted gap depth to merge neighboring peaks (score that can not be exceeded). 
    reverse bool, default False
        Boolean value if False will extend to the right, if True will extend left.

    Returns
    -------
    list :
        Array of extended peak edges.
    """
    if reverse:
        step = -1
    else:
        step = 1

    extended = list()

    peak_num = 0
    for peak in peaks:
        peak_slope = slopes[peak]

        gap = False
        gap_width = 0
        cum_gap_slope = 0

        # start peak extension
        while peak_slope >= min_slope_threshold and peak_slope <= max_slope_threshold or np.abs(gap_width) <= max_gap:
            # stop if extension out of bounds
            if peak + gap_width + step >= len(slopes) or peak + gap_width + step <= 0:
                break

            current_slope = slopes[peak + gap_width]
            current_score = scores[peak + gap_width]
            if not (peak_slope >= min_slope_threshold and peak_slope <= max_slope_threshold) and not gap:
                # start a new gap
                gap = True
                gap_width = -1
                cum_gap_slope = 0
            elif scores[peak] - current_score > gap_depth:
                # end gap and discard; stop peak extension
                break
            elif (cum_gap_slope >= min_slope_threshold and cum_gap_slope <= max_slope_threshold and 
                  current_slope >= min_slope_threshold and current_slope <= max_slope_threshold):
                # end current gap and add it to the peak
                # cum_gap_slope has to be in slope boundaries to ensure equal height between the two plateaus
                # current slope has to be in slope boundaries to ensure new plateau is starting

                gap = False
                peak += gap_width

            # extend
            if gap:
                # extend gap
                gap_width += step
                cum_gap_slope += slopes[peak + gap_width]
            else:
                # extend peak
                peak += step
                peak_slope = slopes[peak]

        extended.append(peak)
        peak_num += 1

    return extended


def merge_peaks(left, right):
    """
    Merge overlapping peaks.

    Parameters
    ----------
    left : list
        List of peak start positions.
    right : list
        List of peak end positions.

    Returns
    -------
    list :
        List of start positions.
    list :
        List of end positions.
    """
    merged_left = list()
    merged_right = list()

    for start, end in zip(left, right):
        if len(merged_left) <= 0:
            merged_left.append(start)

        if len(merged_right) <= 0:
            merged_right.append(end)
        elif merged_right[-1] > start:
            if merged_right[-1] < end:
                merged_right[-1] = end
        else:
            merged_left.append(start)
            merged_right.append(end)

    return merged_left, merged_right


def subtract_footprint(start, end, m_start, m_end):
    """
    Remove m_start, m_end ranges from start, end.

    Parameters
    ----------
    start : list
        List of ranges start positions (beginning of footprint).
    end : list
        List of ranges end positions (end of footprint).
    m_start : list
        List of ranges end positions to remove (beginning of motif sites).
    m_end : list
        List of ranges end positions to remove (end of motif sites).

    Returns
    -------
    list :
        List of start and end positions after substraction.
    """
    subtracted_left = list()
    subtracted_right = list()

    if len(m_start) < 1:
        return [start], [end]

    while len(m_start) > 0:
        cur_m_start = m_start.pop(0)
        cur_m_end = m_end.pop(0)

        if len(subtracted_left) < 1 and start < cur_m_start:
            # add first subtract
            subtracted_left.append(start)
            subtracted_right.append(cur_m_start)

        if len(m_start) > 0:
            # add middle subtracts
            subtracted_left.append(cur_m_end)
            subtracted_right.append(m_start[0])

        if len(m_end) < 1 and cur_m_end < end:
            # add last subtract
            subtracted_left.append(cur_m_end)
            subtracted_right.append(end)

    return subtracted_left, subtracted_right


def plot_stats(peaks, footprints, new_footprints, output='statistics.pdf'):
    """
    Create multiple statistics plots.

    Parameters
    ----------
    peaks : pd.DataFrame
        Dataframe in bed format of the initial peaks (first columns are [chr, start, end]).
    footprints : pd.DataFrame
        Dataframe in bed format of all footprints (needed columns are [chr, start, end, motif_length])
    new_footprints : pd.DataFrame
        Dataframe in bed format of the unknown footprints (first columns are [chr, start, end]).
    output path : str, default statistics.pdf
        Output file.
    """

    fig, axs = plt.subplots(ncols=2, nrows=3, figsize=(14, 15), tight_layout={'rect': (0, 0, 1, 0.95)}) # prevent label clipping; leave space for title
    fig.suptitle('Footprint extraction', fontsize=16)

    # footprint plot
    axs[0][0].set_title('Footprint length distribution')
    footprints["length"] = footprints["end"] - footprints["start"]

    sns.histplot(footprints["length"], ax=axs[0][0], kde=True)

    # new footprint plot
    axs[0][1].set_title('New footprint length distribution')
    new_footprints["length"] = new_footprints["end"] - new_footprints["start"]

    sns.histplot(new_footprints["length"], ax=axs[0][1], kde=True)

    # footprint plot
    axs[1][0].set_title('Percentage coverage of motifs per footprint')
    footprint_coverage = footprints["motif_length"] / footprints["length"]

    sns.histplot(footprint_coverage, ax=axs[1][0], kde=True)

    # nucleotide bar plot
    axs[1][1].set_title('Nucleotide proportions')
    peaks["length"] = peaks.iloc[:, 2] - peaks.iloc[:, 1]

    peak_sum = np.sum(peaks["length"])
    all_sum = np.sum(footprints["length"])
    new_sum = np.sum(new_footprints["length"])

    sns.barplot(x=[peak_sum, all_sum, new_sum], 
                y=["Peaks total nucleotides", "All footprints total nucleotides", "New footprints total nucleotides"],
                ax=axs[1][1])

    # overview count plot
    axs[2][0].set_title('Count')

    sns.barplot(x=[len(peaks), len(footprints), len(new_footprints)],
                y=["Peaks", "Footprints", "New footprints"],
                ax=axs[2][0])

    # footprint binding score distribution plot
    axs[2][1].set_title('Binding score distribution')

    peaks["Origin"] = "Peak"
    footprints["Origin"] = "Footprint"
    new_footprints["Origin"] = "New Footprint"
    score_table = pd.concat([peaks, footprints, new_footprints], ignore_index=True)
    # sns.histplot(data=score_table, x="score", hue="Origin", ax=axs[2][1], kde=True)
    # sns.ecdfplot(data=score_table, x="score", hue="Origin", ax=axs[2][1])
    sns.boxplot(data=score_table, x="score", y="Origin", ax=axs[2][1])

    fig.savefig(output)
    plt.close()


def footprint_extraction(fasta,
                         bigwig,
                         peak_bed,
                         output_dir,
                         motif=None,
                         extend_region=100,
                         min_local_height=0.05,
                         percentile_slope=35,
                         min_size=20,
                         max_size=100,
                         min_score=0.5,
                         gap=10,
                         gap_depth=0.2,
                         moods_threshold=0.00005,
                         log_plot="statistics.pdf",
                         plot_fp=False):
    """
    Extract footprint peaks out of peak regions and filter for footprints of unknown motif origin, keep footprints without a motif site. 

    Parameters
    ----------
    fasta : str
        Fasta file, contains genome for motif binding.
    bigwig : str
        Bigwig file, contains continuouse score for footprint extraction.
    peak_bed : str,
        Bed file, contains peaks from which the footprints are extracted.
    output_dir : str,
        Where the output files should be saved.
    motif : str, default None
        Motif file, contains motifs that the footprints are filtered against.
    extend_region : int, default 100
        Value to extend peak region on both sites.
    min_local_height : float, default 0.05
        Percentage of the maximum peaks height.
    percentile_slope : float, default 35
        Percentage to be added/ subtracted from median (50%) to get upper and lower threshold for footprint extension.
    min_size : int, default 20
        Minimum footprint size.
    max_size : int, default 100
        Maximum footprint size.
    min_score : float, default 0.5
        Minimum footprint score.
    gap : int, default 10
        Maximum gap size between two footprints to be merged.
    gap_depth : float, default 0.2
        Percentage of max footprint height gaps can not exceed.
    moods_threshold : float, default 0.00005
        Sensitivity with which motifs are considered bound (closer to 0 = more sensitive).
    log_plot : str, default statistics.pdf
        File to create statistics summary plots.
    plot_fp : bool, default False
        Boolean if True will create a plot for each peak with footprints. Extremly resource and time intensive and should be avoided for more than 1000 peaks.
    """
    t0 = time.time()

    ##### load and setup #####
    print("Load and setup")
    bigwig = pyBigWig.open(bigwig) # contains scores
    bed = pd.read_csv(peak_bed, delimiter="\t", header=None) # contains peak regions
    bed.sort_values([0,1,2], inplace=True) # probably already sorted but just to be safe

    # extend regions and merge overlapping peaks
    def extend_row(r):
        r[1] = r[1] - extend_region if r[1] - 100 >= 0 else 0
        r[2] = r[2] + extend_region if r[2] + 100 <= bigwig.chroms(r[0]) else bigwig.chroms(r[0])
        return r

    bed = bed.apply(extend_row, axis=1)
    tmp_bed = {'chr': [], 'start': [], 'end': []}
    for chr in list(dict.fromkeys(bed[0])): # make unique chr list while perserving the order
        chr_subset = bed[bed[0] == chr]

        chr_start, chr_end = merge_peaks(chr_subset[1], chr_subset[2])
        tmp_bed['chr'] += [chr] * len(chr_start)
        tmp_bed['start'] += chr_start
        tmp_bed['end'] += chr_end
    bed = pd.DataFrame(tmp_bed)

    # load & setup motifs
    motifs = MotifList().from_file(motif)
    # setup scanner
    for motif in motifs:
        motif.get_threshold(moods_threshold)
        motif.set_prefix()   

    motifs.setup_moods_scanner(".")

    # load genome
    with open(fasta) as fasta:
        fasta_records = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))

    ##### Footprint extraction #####
    plots = list()

    footprint_bed = {'chr': [], 'start': [], 'end': [], 'name': [], 'score': [], 'motif_length': []}
    motif_bed = {'chr': [], 'start': [], 'end': []}
    new_footprint_bed = {'chr': [], 'start': [], 'end': [], 'name': [], 'score': []}

    for index in range(len(bed)):
        ##### find footprints #####
        entry = bed.loc[index]

        chr = entry[0]
        start = entry[1]
        stop = entry[2]

        # extract scores from bigwig
        scores = np.nan_to_num(np.array(bigwig.values(chr, start, stop)))

        # skip peak if all scores are equal
        if len(set(scores)) == 1:
            continue

        # find footprints
        min_prominence = np.max(scores) * min_local_height
        peaks, properties = find_peaks(scores, plateau_size=(min_size, None), prominence=(min_prominence, None), wlen=200)    

        # skip if no footprints found
        if len(peaks) < 1:
            continue

        # interpolate and derive scores
        x = range(0, len(scores))
        spline = UnivariateSpline(x=x, y=scores, s=0)
        slope = spline.derivative()(x)

        ##### extend footprints #####
        # lower threshold
        lt = np.percentile(slope, 50 - percentile_slope)
        # upper threshold
        ut = np.percentile(slope, 50 + percentile_slope)
        # max gap depth
        max_gap_depth = np.max(scores) * gap_depth # score based

        properties["right_edges_extended"] = extend_peak(scores=scores, 
                                                        slopes=slope, 
                                                        peaks=properties["right_edges"].copy(),
                                                        min_slope_threshold=lt,
                                                        max_slope_threshold=ut,
                                                        max_gap=gap,
                                                        gap_depth=max_gap_depth)

        properties["left_edges_extended"] = extend_peak(scores=scores, 
                                                        slopes=slope, 
                                                        peaks=properties["left_edges"].copy(),
                                                        min_slope_threshold=lt,
                                                        max_slope_threshold=ut,
                                                        max_gap=gap,
                                                        gap_depth=max_gap_depth,
                                                        reverse=True)

        # merge overlapping footprints
        properties["left_edges_extended"], properties["right_edges_extended"] = merge_peaks(properties["left_edges_extended"], properties["right_edges_extended"])

        ##### remove motifs #####
        accepted_fp = 0
        accepted_mo = 0
        for f_chr, f_start, f_end in zip([chr] * len(properties["left_edges_extended"]), properties["left_edges_extended"], properties["right_edges_extended"]):
            local_f_start = f_start
            local_f_end = f_end
            f_start += start
            f_end += start

            # filter size and score
            footprint_score = np.mean(scores[local_f_start:local_f_end])
            if footprint_score < min_score or f_end - f_start < min_size or f_end - f_start > max_size:
                continue

            # get sequence
            sequence = str(fasta_records[f_chr][f_start:f_end].seq)

            # scan for motifs, sort and merge hits
            current_hits = motifs.scan_sequence(sequence, OneRegion([f_chr, f_start, f_end]))
            current_hits = sorted(current_hits, key=lambda hit: (hit[1], hit[2]))

            # merge overlapping motif sites
            motif_left, motif_right = merge_peaks([e[1] for e in current_hits], [e[2] for e in current_hits])

            # add to motif dict
            motif_bed['chr'] += [f_chr] * len(motif_left)
            motif_bed['start'].extend(motif_left)
            motif_bed['end'].extend(motif_right)

            # add to footprint dict
            footprint_bed['chr'].append(f_chr)
            footprint_bed['start'].append(f_start)
            footprint_bed['end'].append(f_end)
            footprint_bed['motif_length'].append(np.sum([r - l for l, r in zip(motif_left, motif_right)]))
            footprint_bed['score'].append(footprint_score)
            footprint_bed['name'].append("footprint")

            # remove motif sites from footprints
            new_footprint_left, new_footprint_right = subtract_footprint(start=f_start,
                                                                        end=f_end,
                                                                        m_start=motif_left.copy(),
                                                                        m_end=motif_right.copy())

            # add to new footprint dict
            new_footprint_bed['chr'] += [f_chr] * len(new_footprint_left)
            new_footprint_bed['start'].extend(new_footprint_left)
            new_footprint_bed['end'].extend(new_footprint_right)
            new_footprint_bed['score'].extend([np.mean(scores[l - start:r + 1 - start]) 
                                            for l, r in zip(new_footprint_left, new_footprint_right)])
            new_footprint_bed['name'] += ["new_footprint"] * len(new_footprint_left)

            accepted_fp += 1
            accepted_mo += len(motif_left)

        ##### plot #####
        if plot_fp and accepted_fp > 0:
            # get footprints and motifs of current peak
            # also subtract start to match score index
            footprints_plot = {'start': [i - start for i in footprint_bed['start'][-accepted_fp:]],
                        'end': [i - start for i in footprint_bed['end'][-accepted_fp:]]}
            if accepted_mo > 0:
                motifs_plot = {'start': [i - start for i in motif_bed['start'][-accepted_mo:]],
                            'end': [i - start for i in motif_bed['end'][-accepted_mo:]]}
            else:
                motifs_plot = None

            plot = plot_footprints(scores=scores,
                                slope=slope,
                                percentile_threshold=percentile_slope,
                                footprints=footprints_plot,
                                motifs=motifs_plot,
                                title="    ".join(list(map(str, bed.loc[index])))
                                )

            plots.append(plot)
            plt.close() # prevent plots from showing

    ##### create table and save as bed #####
    footprint_bed = pd.DataFrame(footprint_bed)
    footprint_bed.iloc[:, 0:5].to_csv(os.path.join(output_dir, "footprints.bed"), index=False, header=False, sep="\t")

    motif_bed = pd.DataFrame(motif_bed)
    motif_bed.to_csv(os.path.join(output_dir, "motifs.bed"), index=False, header=False, sep="\t")

    new_footprint_bed = pd.DataFrame(new_footprint_bed)
    # remove footprints below min_size
    new_footprint_bed.drop(new_footprint_bed[new_footprint_bed["end"] - new_footprint_bed["start"] < min_size].index, inplace=True)
    new_footprint_bed.iloc[:, 0:5].to_csv(os.path.join(output_dir, "new_footprints.bed"), index=False, header=False, sep="\t")

    if plot_fp:
        with PdfPages(os.path.join(output_dir, "footprints.pdf")) as pdf:
            for plot in plots:
                pdf.savefig(plot)

    # plot statistics
    if log_plot is not None:
        plot_stats(bed,
                footprint_bed,
                new_footprint_bed,
                output=log_plot)

    print(f"Total runtime {(time.time() - t0) / 60} min.")


if __name__ == "__main__":
    # parse command line args
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Extract footprints by utilizing a binding score and optionally filter for ones not containing a motif site.')

    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument('-f', '--fasta', required=True, help='Input fasta file of the genome.')
    required.add_argument('-b', '--bigwig', required=True, help='Input bigwig file.')
    required.add_argument('-p', '--peak_bed', required=True, help='Input bed file containing the peaks in which footprints will be searched.')
    required.add_argument('-o', '--output_dir', required=True, help='Path to the output directory.')

    optional.add_argument('-m', '--motif', help='Input motif file.', default=None)
    optional.add_argument('-e', '--extend_region', default=100, type=int, help='Extend peaks on both sites by given number of bases.')
    optional.add_argument('-l', '--min_local_height', default=0.05, type=float, help='Percentage of the maximum footprints height.')
    optional.add_argument('-s', '--percentile_slope', default=35, type=int, help='Accepted slope range beginning from 50%% Quantile. E.g. for 35 all slopes between the 15%% and 85%% quantile will be accepted.')
    optional.add_argument('-i', '--min_size', default=20, type=int, help='Minimum size of a footprint.')
    optional.add_argument('-a', '--max_size', default=100, type=int, help='Maximum size of a footprint.')
    optional.add_argument('-r', '--min_score', default=0.5, type=float, help='Minimum score of a footprint. Calculated as the mean of all scores within the footprint.')
    optional.add_argument('-g', '--gap', default=10, type=int, help='Maximum gap width between two footprints to be merged.')
    optional.add_argument('-d', '--gap_depth', default=0.2, type=float, help='Maximum gap depth between to footprints to be merged. Percentage calculated from the maximum footprint height in the peak.')
    optional.add_argument('-c', '--moods_threshold', default=0.00005, type=float, help='P-value like filter for motif binding. Values closer to 0 mean smaller distance between motif and binding site.')
    optional.add_argument('-j', '--log_plot', default="statistics.pdf", help='Create summary plots at given file location. None to disable.')

    args = parser.parse_args()

    # print help if no parameters given
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()

    args = vars(args)

    # run
    footprint_extraction(**args)
