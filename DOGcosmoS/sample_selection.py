import numpy as np
from velociraptor import load as load_catalogue
import configparser

def read_paths_from_config(config_file='specify_your_paths.ini'):
    """Reads the paths from the configuration file 'specify_your_paths.ini'.
    
    If you get an error because the halo_catalogue_files cannot be found, check if your halo catalogue filenames end in '.0' and if not, edit the respective line in this function accordingly.
    
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    
    snapshotfile = config.get('INPUT', 'snapshotfile_path')
    halo_catalogue = config.get('INPUT', 'halo_catalogue_path')
    halo_catalogue_files = dict(
        properties=f"{halo_catalogue}.properties.0",
        catalog_groups=f"{halo_catalogue}.catalog_groups.0",
    )
    output_path = config.get('OUTPUT', 'output_directory_path')

    return snapshotfile, halo_catalogue_files, output_path
    

def select_sample(vmax_min=60, vmax_max=120, centrals=True, save=None):
    """Function to select your galaxy sample according to the range in maximum rotational velocity and environment.

    Args:
        vmax_min: lower limit of maximum rotational velocity that your sample should include. Default = 60 km/s
        vmax_max: upper limit of maximum rotational velocity that you sample should include. Default = 120 km/s.
        centrals: if "True", the sample will be restricted to field galaxies; if "False", the sample will contain satellite galaxies as well. Default = "True"
        save: if "True", the resulting sample indices are saved to your Results directory. Default = "True"
    """
    snapshotfile, catalogue_files, output_path = read_paths_from_config(config_file='specify_your_paths.ini')    
    halos = load_catalogue(catalogue_files['properties'])
    
    if centrals:
        sample_indices = np.argwhere((halos.velocities.vmax > 60) & (halos.velocities.vmax < 120) & (halos.centrals == True)).flatten()
    else:# centrals == "False":
        sample_indices = np.argwhere((halos.velocities.vmax > 60) & (halos.velocities.vmax < 120)).flatten()
        
    if save is not None:
        np.save(output_path+'sample_indices.npy', sample_indices)
    return sample_indices
