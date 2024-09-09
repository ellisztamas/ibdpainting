import numpy as np
import pytest
import ibdpainting as ip

input = 'tests/test_data/panel_to_test.hdf5'
reference = 'tests/test_data/reference_panel.hdf5'
chr1 = 'tests/test_data/reference_panel_chr1.hdf5'

def test_load_genotype_data_gives_right_output():
    """
    Check that load_genotype_data gives the right answers when it should.
    """
    out = ip.load_genotype_data(
        input = input,
        reference = reference,
        sample_name = 'S2.15.002'
    )
    real_names = np.array(['S2.15.002', '1158', '6024', '6184', '8249'])
    assert all(
        [ x == y for x,y in zip(out.samples, real_names) ]
    )

    assert len(out.chr) == 550
    assert len(out.pos) == 550
    assert out.geno.shape == (550, 5, 2)

def test_load_genotype_data_fails_if_missing_sample():
    """
    Check that load_genotype_data fails if the sample name is not in the input file.
    """
    with pytest.raises(Exception):
        ip.load_genotype_data(
            input = input,
            reference = reference,
            sample_name = 'not_a_real_sample_name'
        )

def test_load_genotype_data_fails_if_contigs_dont_match():
    with pytest.raises(Exception):
        ip.load_genotype_data(
            input = input,
            reference = chr1,
            sample_name = '1158_2'
        )
