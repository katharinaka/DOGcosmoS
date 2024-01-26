import numpy as np
from swiftgalaxy import Velociraptor, SWIFTGalaxy
import unyt as U
from velociraptor import load as load_catalogue
from manage_paths import manage_paths
from martini.sources import SWIFTGalaxySource
from martini.sources._L_align import L_align

snapshotfile, halo_catalogue, catalogue_files, path_results, path_plots = manage_paths("COLIBRE_v2023_03")

sample_indices = np.load(path_results+'sample_indices.npy')

def angmom(xyz, vxyz, m, frac):

    transposed = False
    if xyz.ndim != 2:
        raise ValueError(
            "L_align: cannot guess coordinate axis for input with" " ndim != 2."
        )
    elif (xyz.shape[0] == 3) and (xyz.shape[1] == 3):
        raise ValueError(
            "L_align: cannot guess coordinate axis for input with" " shape (3, 3)."
        )
    elif xyz.shape[1] == 3:
        xyz = xyz.T
        vxyz = vxyz.T
        transposed = True
    elif (xyz.shape[0] != 3) and (xyz.shape[1] != 3):
        raise ValueError(
            "L_align: coordinate array shape "
            + str(xyz.shape)
            + " invalid (one dim must be 3)."
        )
        
    rsort = np.argsort(np.sum(np.power(xyz, 2), axis=0), kind="quicksort")
    p = m[np.newaxis] * vxyz
    L = np.cross(xyz, p, axis=0)
    p = p[:, rsort]
    L = L[:, rsort]
    m = m[rsort]
    mcumul = np.cumsum(m) / np.sum(m)
    Nfrac = np.argmin(np.abs(mcumul - frac))
    Nfrac = np.max([Nfrac, 100])  # use a minimum of 100 particles
    Nfrac = np.min([Nfrac, len(m)])  # unless this exceeds particle count
    p = p[:, :Nfrac]
    L = L[:, :Nfrac]
    Ltot = np.sqrt(np.sum(np.power(np.sum(L, axis=1), 2)))
    Lhat = np.sum(L, axis=1) / Ltot
    return Lhat

Lhat_variation_sample =[]
for i in sample_indices:
        
    sg = SWIFTGalaxy(
          snapshotfile,    
             Velociraptor(
             #velociraptor_filebase=halo_catalogue,
             velociraptor_files=catalogue_files,
             halo_index=i,
             centre_type='minpot',
             extra_mask='bound_only',
         ),
    )          
          

    Lhat_varying_massfrac = []

    for frac in [0.1, 0.15, 0.2, 0.25, 0.3]:#, 0.4, 0.5]:#, 0.6, 0.7, 0.8, 0.9]:
          
        Lhat = angmom(sg.gas.coordinates, sg.gas.velocities, sg.gas.atomic_hydrogen_masses, frac)
        Lhat_varying_massfrac.append(Lhat)

    Lhat_variation = np.sqrt(np.mean(np.std(Lhat_varying_massfrac, axis=0)**2))
    #Lhat_variation = np.std(Lhat_varying_massfrac, axis=0)
    Lhat_variation_sample.append(Lhat_variation)
    if Lhat_variation > 0.1:
        print(i, Lhat_variation)

#print(Lhat_variation_sample)
#np.save(path_results+'Lhat_variation', Lhat_variation_sample)
