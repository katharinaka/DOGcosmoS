import sys
sys.path.append('/Users/puk/OSS_development/DOGcosmoS')

import numpy as np
import configparser
from DOGcosmoS.sample_selection import select_sample, read_paths_from_config
from DOGcosmoS.rotation_curve_shapes import quantify_rotation_curve_shapes
from DOGcosmoS.vcirc import vcirc_total, vcirc_particle_type

def test_quantify_rotation_curve_shapes():

    vmax_sample, vfid_sample, eta_rot_sample = quantify_rotation_curve_shapes(vmax_min_sample=50, vmax_max_sample=200, save=True)
    print(vmax_sample)
    print(vfid_sample)

    assert vmax_sample is not None
    assert vfid_sample is not None
    assert eta_rot_sample is not None