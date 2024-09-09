import numpy as np
from warnings import warn
import numpy.ma as ma

class geneticDistance(object):
    """
    A simple class to compare genotype data genetic distances between individuals.  

    Parameters
    ==========
    samples: array
        Vector of length m giving names for each sample.
    chr: array
        Vector of length n giving chromosome labels for each SNP.
    pos: array
        Vector of length n giving SNP positions. Note that SNP positions are inherited from 
        skikit allel and give row numbers from the input HDF5 file rather than
        base-pair positions on each chromosome.
    geno: array
        m x n x 2 array of genotype data where axes index SNPs, individuals, and 
        homologous chromosomes.

    Attributes
    ==========
    samples: array
        Vector of M sample names. The first sample is the input individual to be
        compared to the remaining reference individuals.
    chr: array
        Vector of chromosome labels. These are imported from the reference panel.
    pos: array
        Vector of N SNP positions. These are imported from the reference panel.
    geno: array
        NxMx2 array of genotype data, where N is the number of SNPs and M is the
        number of samples.

    Methods
    =======
    split_into_windows
        Split a geneticDistance object into windows.
    pairwise_distance
        Calculate pairwise genetic distance between an input individual and all 
        reference individuals.
    
    """
    def __init__(self, samples, chr, pos, geno):
        self.samples = samples
        self.chr = chr
        self.pos = pos
        self.geno = geno

    def split_into_windows(self, window_size: int):
        """
        Split a geneticDistance object into windows.

        Splits the geneticDistance object into chromosomes, then into windows on each
        chromosome. It returns a dictionary of geneticDistance objects for each window.

        Parameters
        ==========
        window_size: int
            Window size in base pairs.

        Returns
        =======
        A dictionary of geneticDistance objects with an element for each window.
        Indexes are in the form "Chr:start-stop".
        """
        # Empty dict to store the output
        list_of_distance_objects = {}

        for chr in np.unique(self.chr):
            # Boolean array indexing items in this chromosome
            chr_ix = self.chr == chr
            
            # Array of starting positions for each window. 
            start_positions = np.arange(0, self.pos[chr_ix].max(), window_size)
            for start in start_positions:
                stop  = start + window_size
                # Index positions of SNPs within the current window
                window_ix = (self.pos[chr_ix] >= start) & (self.pos[chr_ix] < stop)
                # Create an object for each window.
                window_name = str(chr) + ":" + str(start) + "-" + str(stop)
                list_of_distance_objects[window_name] = geneticDistance(
                        samples = self.samples,
                        chr  = self.chr[chr_ix][window_ix],
                        pos  = self.pos[chr_ix][window_ix],
                        geno = self.geno[chr_ix][window_ix]
                    )
        
        return list_of_distance_objects
    
    def pairwise_distance(self):
        """
        Calculate pairwise genetic distance between an input individual and all 
        reference individuals.

        The input individual is always the first in the list of samples. Genetic
        distance is the number of allelic differences at each locus between each
        pair, summed over all loci. The calculation is done using masked arrays to
        account for missing data.

        Returns
        =======
        Vector of distances

        """
        masked_geno = ma.masked_array(self.geno, self.geno < 0)

        # Calculate differences at each locus
        per_locus_difference = abs(masked_geno.sum(2)[:,[0]] - masked_geno.sum(2)[:,1:]) / 2
        # Average over loci
        dxy = per_locus_difference.mean(0)
        
        return ma.filled(dxy, -9)
