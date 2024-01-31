import sys
sys.path.append('/Users/puk/OSS_development/DOGcosmoS/DOGcosmoS/')

import numpy as np
import configparser
#from DOGcosmoS import sample_selection
from sample_selection import select_sample, read_paths_from_config 

def test_read_paths():
    snapshotfile, halo_catalogue_file_base, halo_catalogue_files, output_path = read_paths_from_config(config_file='specify_your_paths.ini')

    # insert your paths specified in specify_yout_paths.ini here to assert that they are read correctly.
    assert snapshotfile == '/Users/puk/Downloads/cosmo_volume_example.hdf5'
    assert halo_catalogue_files['properties'] == '/Users/puk/Downloads/cosmo_volume_example.properties'
    assert output_path == '/Users/puk/Test'


def test_select_sample():
    
    sample_indices = select_sample(vmax_min=200, vmax_max=250)

    assert len(sample_indices) > 0
