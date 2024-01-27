import os
import sys
import numpy as np
from swiftgalaxy.halo_finders import Velociraptor
from swiftgalaxy.reader import SWIFTGalaxy
import unyt as u
import astropy.units as U
from velociraptor import load as load_catalogue
from manage_paths import manage_paths
    

def select_sample(vmax_min=60, vmax_max=120, centrals="True", save="True"):

   
    """Function to select your galaxy sample according to the range in maximum rotational velocity and environment.

    Args:
        vmax_min: lower limit of maximum rotational velocity that your sample should include. Default = 60 km/s
        vmax_max: upper limit of maximum rotational velocity that you sample should include. Default = 120 km/s.
        centrals: if "True", the sample will be restricted to field galaxies; if "False", the sample will contain satellite galaxies as well. Default = "True"
        save: if "True", the resulting sample indices are saved to your Results directory. Default = "True"
    """
    
    snapshotfile, halo_catalogue, catalogue_files, path_results, path_plots = manage_paths(simulation_name)
    halos = load_catalogue(catalogue_files['properties'])    
    if centrals == "True":
        sample_indices = np.argwhere((halos.velocities.vmax > 60) & (halos.velocities.vmax < 120) & (halos.centrals == True)).flatten()
    elif centrals == "False":
        sample_indices = np.argwhere((halos.velocities.vmax > 60) & (halos.velocities.vmax < 120)).flatten()
    if save == "True":
        np.save(path_results+'sample_indices.npy', sample_indices)
    return
