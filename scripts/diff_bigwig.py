import pyBigWig
import tqdm
import pandas as pd
import numpy as np

def differential_bigwig(input1, input2, output):
    """
    Compute differential score from the input files. Keeps scores that are higher in input1
    in other words peaks unique or higher in input1.

    Will return input1 score if input1-input2 > threshold, otherwise 0.
    
    Parameters
    ----------
        input1 : str
            Path to bigwig with peaks to filter (minuend).
        input2 : str
            Path to bigwig with peaks to subtract from (subtrahend).
        output : str
            Path to result bigwig.
    """
    # open files
    minuend = pyBigWig.open(input1)
    subtrahend = pyBigWig.open(input2)
    difference = pyBigWig.open(output, "w")
    
    # add file header
    difference.addHeader([(key, val) for key, val in minuend.chroms().items()])
    
    # calculate new scores per chromosome
    chroms = minuend.chroms().keys()
    for chr in tqdm.tqdm(chroms):
        # get minuend and subtrahend info: start, end, score
        info_m = pd.DataFrame(minuend.intervals(chr), columns=["start", "end", "m_val"])
        info_s = pd.DataFrame(subtrahend.intervals(chr), columns=["start", "end", "s_val"])
        
        # stop if minuend is empty
        if not len(info_m):
            continue
        
        # join minuend and subtrahend tables
        info = pd.merge(left=info_m, right=info_s, how="left", on=["start", "end"])
        
        # subtract scores
        info.loc[np.isnan(info["s_val"]), "difference"] = info["m_val"]
        info.loc[np.isnan(info["m_val"]), "difference"] = -info["s_val"]
        info.loc[~(np.isnan(info["m_val"])) & ~(np.isnan(info["s_val"])), "difference"] = info["m_val"] - info["s_val"]
        
        # create new score
        info["score"] = 0.0
        info.loc[(np.isnan(info["difference"])) | (info["difference"] > 0), "score"] = info["m_val"]
        
        # write to output
        difference.addEntries(chr, list(info["start"]), values=list(info["score"]), span=1, step=1)
    
    # close files
    minuend.close()
    subtrahend.close()
    difference.close()