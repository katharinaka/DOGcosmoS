import numpy as np
from swiftgalaxy import Velociraptor, SWIFTGalaxy
import unyt as U
from velociraptor import load as load_catalogue
from manage_paths import manage_paths

simulation_types = ["COLIBRE_v2022_02", "COLIBRE_v2023_03", "COLIBRE_v2023_08", "COLIBRE_v2023_08_lowres", "single-phase", "WDM"]

for simulation_type in simulation_types:
    print(simulation_type)

    result = manage_paths(simulation_type=simulation_type)
    
    # Check if the result is a string (error message)
    if isinstance(result, str):
        print("Error:", result)
    else:
        # Unpack the values only if the result is not an error message
        snapshotfile, halo_catalogue, catalogue_files, path_results, path_plots = result
        sample_indices = np.load(path_results+'sample_indices.npy')
        for i in sample_indices:
            sg = SWIFTGalaxy(
                  snapshotfile,    
                     Velociraptor(
                     velociraptor_files=catalogue_files,
                     halo_index=i,
                     centre_type='minpot',
                     extra_mask='bound_only',
                     ),
                    )
            vgas = sg.gas.velocities
            j