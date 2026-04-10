import h5py
import numpy as np
import textwrap

def require_dataset(f: h5py.File, path: str) -> h5py.Dataset:
    """
    Validate that `path` exists and is a Dataset in an open HDF5 file.
    Returns a lazy h5py.Dataset handle — no data is loaded into memory.
    """
    if path not in f:
        raise KeyError(f"Expected dataset {path!r} not found in {f.filename}.")
    obj = f[path]
    if not isinstance(obj, h5py.Dataset):
        raise TypeError(f"{path!r} is a {type(obj).__name__}, expected an h5py.Dataset.")
    return obj


def _raise_duplicate_marker_error(*, panel_label: str, total: int, unique: int) -> None:
    msg = textwrap.dedent(f"""
        The {panel_label} dataset contains duplicate markers.
        The dataset contains {total} markers, but only {unique} are unique.

        This is usually caused when multiallelic SNPs are collapsed into separate biallelic SNPs.

        Possible culprits:
          - bcftools merge with --merge none
          - Conversion to HDF5 using scikit-allel

        Remove these SNPs and recreate the HDF5 file.
    """).strip()
    raise ValueError(msg)


def find_matching_markers(input_path: str, reference_path: str, sample_name: str) -> dict:
    """
    Prepare HDF5 files for further processing.

    Imports HDF5 files for the test and reference panels and run initial checks.
    This then compares the lists of SNPs in each panel and identifies those that
    are shared.

    Parameters
    ----------
    input_path : str
        Path to a an HDF5 file containing genotype data for one or more samples to check
    reference_path : str
        Path to a HDF5 file containing genotype data for a panel of reference individuals
        to compare the input indivual against.
    sample_name : str
        Sample name for the individual to check. This must be present in the samples
        in the input file.
    
    Returns
    -------
    A dictionary listing:
    samples : np.array
        A vector of strings giving the test individual and all reference
        individuals.
    chr : np.array
        Array of chromosome labels for shared SNPs
    pos : np.array
        Array of positions labels for shared SNPs
    """
    with h5py.File(input_path, "r") as input_hdf5, h5py.File(reference_path, "r") as ref_hdf5:
        input_samples = require_dataset(input_hdf5, "samples").asstr()[...]
        input_chr     = require_dataset(input_hdf5, "variants/CHROM").asstr()[...]
        input_pos     = require_dataset(input_hdf5, "variants/POS")[...]

        #ref_samples   = require_dataset(ref_hdf5, "samples").asstr()[...]
        ref_chr       = require_dataset(ref_hdf5, "variants/CHROM").asstr()[...]
        ref_pos       = require_dataset(ref_hdf5, "variants/POS")[...]
        if sample_name not in input_samples:
            raise ValueError(f"Sample {sample_name!r} not found in input file {input_path!r}.")

        # Find the position of the individual to test
        sample_ix = int(list(input_samples).index(sample_name))

        # Check that contig labels match (set-wise + ordering)
        input_contigs = np.unique(input_chr)
        ref_contigs = np.unique(ref_chr)
        if input_contigs.shape != ref_contigs.shape or np.any(input_contigs != ref_contigs):
            raise ValueError("Contig labels do not match between the input and reference files.")

        # Concatenate chromosome labels and SNP positions
        snp_names_input = [f"{c}:{p}" for c, p in zip(input_chr, input_pos)]
        snp_names_ref   = [f"{c}:{p}" for c, p in zip(ref_chr, ref_pos)]

        # Check for duplicate SNP positions
        n_in = len(snp_names_input)
        n_in_unique = len(set(snp_names_input))
        if n_in != n_in_unique:
            _raise_duplicate_marker_error(panel_label="input", total=n_in, unique=n_in_unique)

        n_ref = len(snp_names_ref)
        n_ref_unique = len(set(snp_names_ref))
        if n_ref != n_ref_unique:
            _raise_duplicate_marker_error(panel_label="reference", total=n_ref, unique=n_ref_unique)

        # Find markers common to both datasets
        shared = set(snp_names_input) & set(snp_names_ref)

        keep_input = np.fromiter((x in shared for x in snp_names_input), dtype=bool, count=n_in)
        keep_ref   = np.fromiter((x in shared for x in snp_names_ref), dtype=bool, count=n_ref)

    return {"sample_ix": sample_ix, "input": keep_input, "ref": keep_ref}
