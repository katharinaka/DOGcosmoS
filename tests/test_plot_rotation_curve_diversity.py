import sys
sys.path.append('/Users/puk/OSS_development/DOGcosmoS')

import numpy as np
import configparser
import matplotlib
from DOGcosmoS.sample_selection import select_sample, read_paths_from_config
from DOGcosmoS.rotation_curve_shapes import quantify_rotation_curve_shapes
from DOGcosmoS.vcirc import vcirc_total, vcirc_particle_type
from DOGcosmoS.plots.plot_rotation_curve_diversity import plot_vfid_vmax

def test_plot_vfid_vmax():

    plot_vfid_vmax() 

    assert True