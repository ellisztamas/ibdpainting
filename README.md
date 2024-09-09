# ibdpainting

`ibdpainting` is a Python tool to visually validate the identity of crossed individuals
from genetic data.

## Contents

- [Premise](#premise)
- [Installation](#installation)
- [Input data files](#input-data-files)
- [Usage](#usage)
- [Author information](#author-information)
- [Contributing](#contributing)

## Premise

`ibdpainting` addresses the situation where you have multiple individuals
derived from a crosses between individuals in a reference panel, and you want to
verify that the crosses really are the genotype you think they are. Taking the
simple example of a biparental cross, you would expect an offspring of the F2 
generation or later to be a mosaic of regions homozygous for either parent, 
potentially interspersed with heterozygous regions, depending on the generation.
`ibdpainting` is a tool to visualise this mosaic pattern.

## Installation

Install with `pip`:
```
pip install ibdpainting
```

## Input data files

The program requires two HDF5 files created from VCF files:

* **Input panel**: An HDF5 file containing SNPs for the crossed individual(s).
This can contain multiple individuals, but the program will only work on one at
a time.
* **Reference panel**: An HDF5 file conataining SNP information for a panel of reference candidate
parents.

The reason for using HDF5 is that it allows for loading data in chunks,
which is much quicker than loading an entire VCF file into memory every time you
want to check a single sample. I recommend creating this using
[vcf_to_hdf5](https://scikit-allel.readthedocs.io/en/latest/io.html#allel.vcf_to_hdf5)
from `scikit-allel`. For example:
```
import allel
allel.vcf_to_hdf5('example.vcf', 'example.h5', fields='*', overwrite=True)
```

Tips for preparing the data:

* `ibdpainting` will only compare SNPs that intersect the input and reference files.
One one hand, this means that it does not matter if the offspring and reference
files contain SNPs that do not match exactly.
On the other, this may cause problems if you are comparing samples with *loads*
of structural variation.
* It is better to have a smaller number of reliable SNPs than a larger number of 
dubious SNPs. For example, in *Arabopidopsis thaliana* that means only using 
common SNPs located in genes.
* `ibdpainting` creates a subplot for every contig label in the input/reference
panel. If you work on an organism with many chromosomes or incompletely assembled
contigs, this could get messy. There is currently no way to subset which 
contigs are shown, so it is probably easiest to supply input data based on only 
a subset of contigs. The longest contigs are likely to be most informative
because you are more likely to be able to spot recombination break points.

## Usage

After installing, `ibdpainting` can be run as a command line tool as follows

```
ibdpainting \
    --input input_file.hdf5 \
    --reference reference_panel.hdf5 \
    --window_size 500000 \
    --sample_name "my_cross" \
    --expected_match "mother" "father" \
    --outdir path/to/output/directory
```

Explanation of the parameters (see also `ibdpainting --help`):

* `--input`: HDF5 file containing the crossed individuals. See [above](#input-data-files).
* `--reference`: HDF5 file containing the reference panel. See [above](#input-data-files).
* `--window_size`: Window size in base pairs.
* `--sample_name`: Name of the crossed individual to compare to the reference 
panel. This must be present in the input file - you can check the original VCF file with something
like `bcftools query -l $input_vcf.vcf.gz | grep "my_cross"`.
* `--expected_match`: List of one or more expected parents of the test individual.
These names should be among the samples in the reference panel. Names should be
separated by spaces.
* `--outdir`: Path to the directory to save the output.

Additional optional parameters:

* `--keep_ibd_table`: Write an intermediate text file giving genetic distance 
between the crossed individual and each candidate at each window in the genome.
Defaults to False, because these can be quite large.
* `--max_to_plot`: Integer number of candidates to plot.
`ibdpainting` makes an intial ranking of candidates based on genome-wide 
similarity to the test individual, and plots only the top candidates.

## Author information

Tom Ellis

## Contributing

I will repeat the following from the [documentation](https://scikit-allel.readthedocs.io/en/stable/) for `scikit-allel`:

> This is academic software, written in the cracks of free time between other commitments, by people who are often learning as we code. We greatly appreciate bug reports, pull requests, and any other feedback or advice. If you do find a bug, we’ll do our best to fix it, but apologies in advance if we are not able to respond quickly. If you are doing any serious work with this package, please do not expect everything to work perfectly first time or be 100% correct. Treat everything with a healthy dose of suspicion, and don’t be afraid to dive into the source code if you have to. Pull requests are always welcome.
