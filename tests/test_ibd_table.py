import ibdpainting as ip
import numpy as np

input = 'tests/test_data/panel_to_test.hdf5'
reference = 'tests/test_data/reference_panel.hdf5'
ref_vcf = 'tests/test_data/reference_panel.hdf5'
chr1 = 'tests/test_data/reference_panel_chr1.hdf5'

def test_ibd_table():
    ibd = ip.ibd_table(
        input=ref_vcf,
        reference=reference,
        sample_name='1158',
        window_size=1000
    )
    # Check the dataframe is the right shape
    assert ibd.shape == (200, 5)
    # Check the column names are as expected
    assert ibd.keys().to_list() == ['window', '1158', '6024', '6184', '8249']
    # Check that the column for the true parent is all zeroes or -9
    assert all(
        (ibd['1158'] == 1) | (np.isnan(ibd['1158']))
    )
    # Check that a non-parent are not all -9.
    assert any(ibd['8249'] != 0)


def test_ibdtable_returns_expected_F1():
    ibd = ip.ibd_table(
        input=input,
        reference=reference,
        sample_name='F2.05.015',
        expected_match=['6024', '8249'],
        window_size=1000
    )

    # Check the column names are as expected
    assert ibd.keys().to_list() == ['window', '1158', '6024', '6184', '8249', 'Expected_F1']
