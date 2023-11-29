import numpy as np
from swiftgalaxy import Velociraptor, SWIFTGalaxy
import unyt as u
import astropy.units as U
from velociraptor import load as load_catalogue
from manage_paths import manage_paths

simulation_types = ["COLIBRE_v2022_02", "COLIBRE_v2023_03", "COLIBRE_v2023_08", "COLIBRE_v2023_08_lowres", "single-phase", "WDM", "ZOOM_Jemima"]

for simulation_type in simulation_types[6:]:
    print(simulation_type)

    result = manage_paths(simulation_type=simulation_type)
    
    # Check if the result is a string (error message)
    if isinstance(result, str):
        print("Error:", result)
    else:
        # Unpack the values only if the result is not an error message
        snapshotfile, halo_catalogue, catalogue_files, path_results, path_plots = result

    halos = load_catalogue(catalogue_files['properties'])
    
    if simulation_type == "ZOOM_Jemima":
        center = [101.0939839,  182.84324533, 370.79676284]*u.Mpc
        positions = np.concatenate([halos.positions.xcminpot[:,np.newaxis]-center[0], halos.positions.ycminpot[:,np.newaxis]-center[1], halos.positions.zcminpot[:,np.newaxis]-center[2]], axis=1)#.to('Mpc')
        
        r = np.sqrt(np.sum(np.power(positions, 2), axis=1))#.to('Mpc')
        print(r)
        mask_highres_region = np.where(r <= 5*u.Mpc)
        sample_indices = np.argwhere((halos.velocities.vmax[mask_highres_region] > 60) 
                                     & (halos.velocities.vmax[mask_highres_region] < 120) 
                                     & (halos.centrals[mask_highres_region] == True)).flatten()
    else:
        sample_indices = np.argwhere((halos.velocities.vmax > 60) & (halos.velocities.vmax < 120) & (halos.centrals == True)).flatten()
        
    print(sample_indices)
    print(path_results)
    

    np.save(path_results+'sample_indices.npy', sample_indices)
