# Not run - need to work out how to test a CLI properly.

import ibdpainting as ip


input = 'tests/test_data/panel_to_test.hdf5'
reference = 'tests/test_data/reference_panel.hdf5'
ref_vcf = 'tests/test_data/reference_panel.hdf5'
chr1 = 'tests/test_data/reference_panel_chr1.hdf5'


"""
ibdpainting \
    --input tests/test_data/reference_panel.hdf5 \
    --reference tests/test_data/reference_panel.hdf5 \
    --sample_name 1158 \
    --window_size 1000 \
    --outdir tests/test_output \
    --expected_match 1158 \
    --no-interactive
"""