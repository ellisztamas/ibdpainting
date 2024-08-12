import numpy as np
import pandas as pd

def ibd_scores(ibd_table):
    """
    Mean IBD across the genome.

    Calculate mean genetic distance from a test individual to each of a panel of
    reference samples, ignoring windows where there was only missing data.

    Parameters
    ==========
    ibd_table: pd.DataFrame
        DataFrame with a row for each window in the genome and a column for each 
        sample in the reference panel. Elements show genetic distance between the 
        test individual and each reference individual in a single window.
        This is generated by ibd_table().

    Returns
    =======
    A DataFrame with a row for each candidate in the reference panel, and a 
    column indicating mean genetic distance over windows across the genome.
    Values closer to zero indicate that the sample is more likely to be a match.
    """
    # Coerce missing data to NaN for correct column means.
    ibd_table = ibd_table.replace(-9,np.NaN)

    # Get column-mean IBD for each candidate, allowing for missing data
    ibd_scores_for_each_candidate = np.array(
        [ np.nanmean(ibd_table[col]) for col in ibd_table.keys()[1:] ]
    )
    scores = pd.DataFrame({
        'candidate': ibd_table.keys()[1:],
        'score' : ibd_scores_for_each_candidate
    })
    return scores