# Sample files for running unit tests

## VCF files

Sample VCF files for comparing genotyped samples against a reference panel:
```
reference_panel.vcf.gz
panel_to_test.vcf.gz
reference_panel_chr1.vcf.gz
```
These are taken from the methylation crosses project in
`/groups/nordborg/projects/epiclines/002.pedigree.
The reference panel is the four parents of the crosses performed (1158, 6024, 
6184, 8249) taken from Fernandos VCF.
The test panel are 10 samples of parents, F1s, F2s, F3s, a negative control and 
one positive control(Columbia).

The first two files should be the first 100kb on chromosomes 1 and 2.
The third is the same, but chromosome 1 only, to check that the program fails
if chromsome labels don't match.

Code to create the files:
```
bcftools view \
    --regions Chr1:1-100000,Chr2:1-100000 \
    --samples 1158,6024,6184,8249 \
    --output tests/test_data/reference_panel.vcf.gz \
    /groups/nordborg/projects/epiclines/002.pedigree/03_processing/05_validate_genotyping/output/vcf/parents_only.vcf.gz

bcftools view \
    --regions Chr1:1-100000,Chr2:1-100000 \
    --samples F2.05.015,S2.06.002,S2.15.002,F2.27.020 \
    --output tests/test_data/panel_to_test.vcf.gz \
    /scratch-cbe/users/thomas.ellis/meth_pedigree/05_validate_genotyping/02_genotype_calls_by_plate/2021-014.vcf.gz

bcftools view \
    --regions Chr1:1-100000 \
    --samples 1158,6024,6184,8249 \
    --output tests/test_data/reference_panel_chr1.vcf.gz \
    /groups/nordborg/projects/epiclines/002.pedigree/03_processing/05_validate_genotyping/output/vcf/parents_only.vcf.gz

```

## HDF5 files

HDF5 file matching the VCF files described above.

```
import allel

allel.vcf_to_hdf5(
    'tests/test_data/reference_panel.vcf.gz',
    'tests/test_data/reference_panel.hdf5',
    fields='*', overwrite=True
    )

allel.vcf_to_hdf5(
    'tests/test_data/panel_to_test.vcf.gz',
    'tests/test_data/panel_to_test.hdf5',
    fields='*', overwrite=True
    )

allel.vcf_to_hdf5(
    'tests/test_data/reference_panel_chr1.vcf.gz',
    'tests/test_data/reference_panel_chr1.hdf5',
    fields='*', overwrite=True
    )

```

