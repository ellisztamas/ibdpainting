# ibdpainting

`ibdpainting` is a Python tool to visually validate the identity of crossed individuals
from genetic data.

## Contents

- [Premise](#premise)
- [Installation](#installation)
- [Input data files](#input-data-files)
- [Usage](#usage)
- [Output and interpretation](#output-and-interpretation)
- [Author information](#author-information)
- [Contributing](#contributing)

## Premise

`ibdpainting` addresses the situation where you have multiple individuals
derived from a crosses between individuals in a reference panel, and you want to
verify that the crosses really are the genotype you think they are. Taking the
simple example of a biparental cross, you would expect an offspring of the F2 
generation or later to be a mosaic of regions that are identical by descent (IBD)
to either parent, potentially interspersed with heterozygous regions, depending
on the generation.
`ibdpainting` is a tool to visualise this mosaic pattern by 'painting' the
pattern of IBD across the genome.

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

Explanation of the parameters:

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

See the output of `ibdpainting --help` for additional optional arguments.

## Output and interpretation

By default, `ibdpainting` creates three files:

* A `plot_ibd.png` image of the genome, showing the position along each chromomosome
along the x-axis, and the genetic distance from the progeny to each candidate 
along the y-axis. If a candidate parent is IBD to the progeny, points on the 
y-axis should be zero, genotyping errors notwithstanding. Candidate parents
given as expected parents will be shown with coloured lines. The (usually ten)
next-closest other candidates are shown in grey.
* An `ibd_score.csv` file listing possible combinations of candidate parents
and a score for each. The score for a single pair is calculated by the minimum
distance between the offspring and either candidate in each window, and averaging
these over all non-NA windows. Scores close to zero indicate a better match.
A good match will ideally be an order of magnitude better than the next pair.
Only scores for pairs with the 100 most likely candidates are shown.
* `plot_ibd.html`: An optional interactive version of the `png` file. Roll over points to see which candidate is which. These files are about ten times larger than the `png` files. Disable with `--no-interactive`.
* `ibd_table.csv`: An optional text file giving genetic distances from the progeny individual to every candidate in every window. These files are typically big, so are not created by default. Enable with `--keep_ibd_table`.

For examples of the output and the interpretation of different patterns, see [the examples here](https://github.com/ellisztamas/ibdpainting/blob/main/setup.py/example_results.html).

## Author information

Tom Ellis

## Contributing

I will repeat the following from the [documentation](https://scikit-allel.readthedocs.io/en/stable/) for `scikit-allel`:

> This is academic software, written in the cracks of free time between other commitments, by people who are often learning as we code. We greatly appreciate bug reports, pull requests, and any other feedback or advice. If you do find a bug, we’ll do our best to fix it, but apologies in advance if we are not able to respond quickly. If you are doing any serious work with this package, please do not expect everything to work perfectly first time or be 100% correct. Treat everything with a healthy dose of suspicion, and don’t be afraid to dive into the source code if you have to. Pull requests are always welcome.
