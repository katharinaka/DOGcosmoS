import sys
sys.path.append('/Users/puk/OSS_development/DOGcosmoS/DOGcosmoS')

import numpy as np
import configparser
from select_sample import sample_selection, read_paths_from_config 

def test_read_paths():
    snapshotfile, halo_catalogue_files, output_path = read_paths_from_config('DOGcosmoS/specify_your_paths.ini')


    assert snapshotfile == '/Users/puk/Downloads/small_cosmo_volume/snap_0199.hdf5'
    assert halo_catalogue_files is not None
    assert output_path is not None


def test_select_sample():
    
    sample_indices = select_sample(vmax_min=200, vmax_max=250)

    assert len(sample_indices) > 0
