import pandas as pd

from ibdpainting import load_genotype_data

def ibd_table(input:str, reference:str, sample_name:str, window_size:int):
    """
    Compare allele sharing across the genome.

    Calculate genetic distance between a test individual and a panel of
    reference genomes.

    Parameters
    ==========
    input: str
        Path to an HDF5 file containing genotype data for one or more samples to
        test
    reference: str
        Path to an HDF5 file containing genotype data for a panel of reference
        individuals to compare the test individual against.
    sample_name: str
        Sample name for the individual to check.
        This must be present in the samples in the input HDF5 file.
    window_size: int
        Window size in base pairs.

    Returns
    =======
    DataFrame with a row for each window in the genome and a column for each 
    sample in the reference panel. Elements show genetic distance between the 
    test individual and each reference individual in a single window.
    """
    genetic_distance = load_genotype_data(
            input = input,
            reference = reference,
            sample_name = sample_name
        )
    # Divide the genome into windows
    distances_in_windows = genetic_distance.split_into_windows(window_size)

    # Dataframe with a row for each window across the genome and a column for each sample in the reference panel.
    distance_array = pd.DataFrame(
        [ v.pairwise_distance() for v in distances_in_windows.values() ],
        columns = genetic_distance.samples[1:]
    )
    distance_array.insert(0, 'window', distances_in_windows.keys())

    return distance_array