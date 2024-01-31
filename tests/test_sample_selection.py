import sys
sys.path.append('/Users/puk/OSS_development/DOGcosmoS/DOGcosmoS/')

import numpy as np
import configparser
#from DOGcosmoS import sample_selection
from sample_selection import select_sample, read_paths_from_config 

def test_read_paths():
    snapshotfile, halo_catalogue_file_base, halo_catalogue_files, output_path = read_paths_from_config(config_file='specify_your_paths.ini')


    assert snapshotfile == '/Users/puk/Downloads/cosmo_volume_example.hdf5'
    assert halo_catalogue_files['properties'] == '/Users/puk/Downloads/cosmo_volume_example.properties'
    assert output_path is not None


def test_select_sample():
    
    sample_indices = select_sample(vmax_min=200, vmax_max=250)

    assert len(sample_indices) > 0
