import ibdpainting as ip
import numpy as np

hdf5_path = 'tests/test_data/reference_panel.hdf5'

def test_find_matching_markers():
    markers = ip.find_matching_markers(
        input_path=hdf5_path,
        reference_path=hdf5_path,
        sample_name='1158'
    )
    assert isinstance(markers, dict)
    assert list(markers.keys()) == ['sample_ix', 'input', 'ref']
    assert markers['sample_ix'] == 0
    # Check that all markers are valid in this case
    assert all(markers['input'])
    assert all(markers['ref'])
