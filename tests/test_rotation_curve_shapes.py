import sys
sys.path.append('/Users/puk/OSS_development/DOGcosmoS/DOGcosmoS/')

import numpy as np
import configparser
from sample_selection import select_sample, read_paths_from_config
from rotation_curve_shapes import quantify_rotation_curve_shapes
from vcirc import vcirc_total, vcirc_particle_type

def test_quantify_rotation_curve_shapes():

    vmax_sample, vfid_sample, eta_rot_sample = quantify_rotation_curve_shapes(vmax_min_sample=200, vmax_max_sample=250, save=True)
    print(vmax_sample)
    print(vfid_sample)

    assert len(vmax_sample) > 0
    assert len(vfid_sample) > 0
    assert len(eta_rot_sample) > 0