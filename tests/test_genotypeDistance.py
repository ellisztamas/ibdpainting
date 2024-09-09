import ibdpainting as ip


input = 'tests/test_data/panel_to_test.hdf5'
reference = 'tests/test_data/reference_panel.hdf5'
ref_vcf = 'tests/test_data/reference_panel.hdf5'
chr1 = 'tests/test_data/reference_panel_chr1.hdf5'


def test_split_into_windows_functions():
    vcfd = ip.load_genotype_data(
            input = input,
            reference = reference,
            sample_name = 'S2.15.002'
        )
    split_vcfd = vcfd.split_into_windows(1000)
    assert all( split_vcfd['Chr1:0-1000'].pos >= 0 )
    assert all( split_vcfd['Chr1:0-1000'].pos < 1000 )
    assert all(split_vcfd['Chr1:0-1000'].chr == "Chr1")
    assert len(split_vcfd['Chr1:0-1000'].geno.shape) == 3
    # Check you get only one window per chr if window size >> chr length
    assert len(vcfd.split_into_windows(1000000)) == 2

def test_pairwise_distance_works():
    """
    There are four accessions in the reference panel.
    Test each against the whole panel, and check that one of them comes out as
    identical in each case.
    """
    # 1158
    check_1158 = ip.load_genotype_data(
        input = ref_vcf,
        reference = reference,
        sample_name= '1158'
        ).pairwise_distance()

    assert check_1158[0] == 0
    assert all(check_1158[1:] > 0)

    # 6024
    check_6024 = ip.load_genotype_data(input = ref_vcf, reference = reference,
                sample_name= '6024'
        ).pairwise_distance()

    assert check_6024[1] == 0
    assert all(check_6024[[0,2,3]] > 0)

    # 6184
    check_6184 = ip.load_genotype_data(input = ref_vcf, reference = reference,
                sample_name= '6184'
        ).pairwise_distance()

    assert check_6184[2] == 0
    assert all(check_6184[[0,1,3]] > 0)

    # 8249
    check_8249 = ip.load_genotype_data(input = ref_vcf, reference = reference,
                sample_name= '8249'
        ).pairwise_distance()

    assert check_8249[3] == 0
    assert all(check_8249[:2] > 0)
    
    
def test_missing_data_in_geneticDistance():
    """Test that pairwise_distance returns -9 if all loci are NA.
    """
    vcfd = ip.load_genotype_data(
        input = input,
        reference = reference,
        sample_name = 'S2.15.002'
    )   
    vcfd.geno[:,1] = -1
    assert vcfd.pairwise_distance()[0] == -9