import allel
import h5py
import numpy as np
import ibdpainting.geneticDistance

def load_genotype_data(input, reference, sample_name):
    """
    Import and merge test and reference data files.

    Import genotype data for one or more input samples and a panel of reference samples
    to compare to. Subset each so that the markers are really identical. Merge the
    arrays of genotype calls so that the data for the input appear first on the
    first axis of the genotype call arrays.

    Parameters
    ==========
    input: str
        Path to a an HDF5 file containing genotype data for one or more samples to check
    reference: str
        Path to a HDF5 file containing genotype data for a panel of reference individuals
        to compare the input indivual against.
    sample_name: str
        Sample name for the individual to check. This must be present in the samples
        in the input file.

    Return
    ======
    An object of class geneticDistance.
    """
    # Read in the data files
    input_hdf5 = h5py.File(input, mode='r')
    ref_hdf5  = h5py.File(reference, mode="r")

    input_str_data = {
        'samples' : [ x.decode('utf-8') for x in input_hdf5['samples'][:] ],
        'chr'     : [ x.decode('utf-8') for x in input_hdf5['variants/CHROM'][:] ]
    }

    ref_str_data = {
        'samples' : [ x.decode('utf-8') for x in ref_hdf5['samples'][:] ],
        'chr'     : [ x.decode('utf-8') for x in ref_hdf5['variants/CHROM'][:] ]
    }
    import time

    if sample_name not in input_str_data['samples']:
        raise ValueError("The sample name is not in the list of samples in the input file.")
    else: 
        # Find the position of the individual to test
        sample_ix = np.where(
            [ sample_name == x for x in input_str_data['samples'] ]
            )[0][0]
        # Join vectors of sample names, with the test individual first
        new_samples = np.append(
            input_str_data['samples'][sample_ix], ref_str_data['samples']
            )        

    # Check that contig labels match
    chr_labels = {
        'input' : np.unique(input_str_data['chr']),
        'ref'   : np.unique(ref_str_data['chr'])
    }
    if len(chr_labels['input']) != len(chr_labels['ref']):
        raise ValueError(
            "The number of unique contig labels do not match: the input an HDF5 has {}, but the reference panel has {}.".
            format( chr_labels['input'], chr_labels['ref'] )
        )
    elif any( chr_labels['input'] != chr_labels['ref'] ):
        raise ValueError(
            "Contig labels do not match between the input and reference files."
        )
    
    # Make sure we only compare SNPs that are found in both datasets.
    # Concatenate chromosome labels and SNP positions
    snp_names = {
        'input' : [ str(chr) + ":" + str(pos) for chr,pos in zip(input_str_data['chr'], input_hdf5['variants/POS'][:]) ],
        'ref'   : [ str(chr) + ":" + str(pos) for chr,pos in zip(ref_str_data['chr'], ref_hdf5['variants/POS'][:]) ]
    }
    # Find the SNP position names that are common to both datasets
    matching_SNPs_in_both_files = np.intersect1d(
        snp_names['input'],
        snp_names['ref']
        )
    which_SNPs_to_keep = {
        "input" : [ x in matching_SNPs_in_both_files for x in snp_names['input'] ],
        "ref"   : [ x in matching_SNPs_in_both_files for x in snp_names['ref'] ]
    }


    # Append the genotype data for the test individual to the array of the reference panel
    new_geno = np.concatenate(
        (input_hdf5['calldata/GT'][which_SNPs_to_keep['input'], sample_ix][:, np.newaxis],
        ref_hdf5['calldata/GT'][which_SNPs_to_keep['ref']]),
        axis=1
        )
    
    # Define an output before closing the Hdf5 file
    output = ibdpainting.geneticDistance(
        samples = new_samples,
        chr = np.array(ref_str_data['chr'])[np.where(which_SNPs_to_keep['ref'])[0]],
        pos = ref_hdf5['variants/POS'][:][which_SNPs_to_keep['ref']],
        geno = new_geno
    )

    ref_hdf5.close()

    return output