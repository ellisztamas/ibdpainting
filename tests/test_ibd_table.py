import ibdpainting as ip


input = 'tests/test_data/panel_to_test.vcf.gz'
reference = 'tests/test_data/reference_panel.hdf5'
ref_vcf = 'tests/test_data/reference_panel.vcf.gz'
chr1 = 'tests/test_data/reference_panel_chr1.vcf.gz'

def test_ibd_table():
    ibd = ip.ibd_table(
        input=ref_vcf,
        reference=reference,
        sample_name='1158',
        window_size=1000
    )
    # Check the dataframe is the right shape
    assert ibd.shape == (200, 5)
    # Check that the column for the true parent is all zeroes or -9
    assert all(
        (ibd['1158'] == 0) | (ibd['1158'] == -9)
    )
    # Check that a non-parent are not all -9.
    assert any(ibd['8249'] != 0)