import sys
sys.path.append('/Users/puk/OSS_development/DOGcosmoS/DOGcosmoS')

import numpy as np
import configparser
import matplotlib
from sample_selection import select_sample, read_paths_from_config
from rotation_curve_shapes import quantify_rotation_curve_shapes
from vcirc import vcirc_total, vcirc_particle_type
from plot_rotation_curve_diversity import plot_vfid_vmax

def test_plot_vfid_vmax():

    plot_vfid_vmax(200, 250) 

    assert True