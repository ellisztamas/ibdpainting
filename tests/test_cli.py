# Not run - need to work out how to test a CLI properly.

import ibdpainting as ip


input = 'tests/test_data/panel_to_test.vcf.gz'
reference = 'tests/test_data/reference_panel.hdf5'
ref_vcf = 'tests/test_data/reference_panel.vcf.gz'
chr1 = 'tests/test_data/reference_panel_chr1.vcf.gz'


"""
ibdpainting \
    --input tests/test_data/reference_panel.vcf.gz \
    --reference tests/test_data/reference_panel.vcf.gz \
    --sample_name 1158 \
    --window_size 1000 \
    --outdir tests/test_output \
    --expected_match 1158 \
    --no-interactive
"""