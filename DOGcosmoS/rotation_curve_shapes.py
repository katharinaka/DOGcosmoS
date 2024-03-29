import numpy as np
import configparser
from scipy.interpolate import splrep, splev
from sample_selection import read_paths_from_config, select_sample
from vcirc import vcirc_total

def quantify_rotation_curve_shapes(vmax_min_sample, vmax_max_sample, save=None):
#def quantify_rotation_curve_shapes(sample_indices, save=None):
    """Calculates toration curve shape parameters for your galaxy sample.
    
    Arguments:
             vmax_min_sample: Lower limit of maximum rotation velocity that your sample should contain.
             vmax_max_sample: Upper limit of maximum rotation velocity that your sample should contain.
             save:            If True, an numpy array containing vmax and vfid of your sample is saved to your output path. Default is None.
    Returns:
        vmax_sample, vfid_sample, eta_rot_sample
    """
    sample_indices = select_sample(vmax_min=vmax_min_sample, vmax_max=vmax_max_sample, centrals=True, save=None)
    vmax_sample = np.zeros(len(sample_indices))
    rfid_sample = np.zeros(len(sample_indices))
    vfid_sample = np.zeros(len(sample_indices))
    eta_rot_sample = np.zeros(len(sample_indices))
    
    for i in range(len(sample_indices)):
        
        halo_index = sample_indices[i]
        
        #if halo_index == 2:
        #    continue
                
        r, vcirc = vcirc_total(halo_index, save=None, plot=None)
        if not len(vcirc) > 1:
            continue
    
        vmax = np.max(vcirc)
        rfid = 2 * vmax/70   #/70km/s kpc

        # Perform spline interpolation
        spl = splrep(r, vcirc)

        # Evaluate the spline at desired points
        vfid = splev(rfid, spl)
    
        vmax_sample[i] = vmax
        rfid_sample[i] = rfid
        vfid_sample[i] = vfid
        eta_rot_sample[i] = vfid/vmax

    vmax_vfid = np.concatenate([vmax_sample[:,np.newaxis], vfid_sample[:,np.newaxis]], axis=1)
    print('vfid_vmax: ', vmax_vfid)
    
    if save is not None:
        snapshotfile, halo_catalogue_filebase, halo_catalogue_properties, output_path = read_paths_from_config(config_file='specify_your_paths.ini')
        np.save(output_path+'vmax_vfid.npy', vmax_vfid)
        
    return vmax_sample, vfid_sample, eta_rot_sample
