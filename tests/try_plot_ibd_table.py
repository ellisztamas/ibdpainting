import ibdpainting as ip

input = 'tests/test_data/panel_to_test.hdf5'
reference = 'tests/test_data/reference_panel.hdf5'
ref_vcf = 'tests/test_data/reference_panel.hdf5'
chr1 = 'tests/test_data/reference_panel_chr1.hdf5'

ibd = ip.ibd_table(
    input=ref_vcf,
    reference=reference,
    sample_name='1158',
    expected_match = ['1158', '8249'],
    window_size=1000
)

fig = ip.plot_ibd_table(
    ibd_table = ibd,
    sample_name = '1158',
    expected_match = ['1158'],
    max_to_plot=10
)
fig.show()

fig = ip.plot_ibd_table(
    ibd_table = ibd,
    sample_name = '1158',
    expected_match = ['1158'],
    max_to_plot=10,
    plot_heterozygosity = True
)
fig.show()

fig = ip.plot_ibd_table(
    ibd_table = ibd,
    sample_name = '1158',
    expected_match = ['1158', '8249'],
    max_to_plot=10
    )
fig.show()