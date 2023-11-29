import numpy as np
from swiftgalaxy import Velociraptor, SWIFTGalaxy
import unyt as U
from velociraptor import load as load_catalogue
from manage_paths import manage_paths
import matplotlib.pyplot as plt

simulation_types = ["COLIBRE_v2022_02", "COLIBRE_v2023_03", "COLIBRE_v2023_08", "COLIBRE_v2023_08_lowres", "single-phase", "WDM", "ZOOM_Jemima"]

def HI_disk_size(xyz, m_HI, frac=0.9, save=None):
    transposed = False
    if xyz.ndim != 2:
        raise ValueError(
            "HI_disk_size: cannot calculate disk size for input with" " ndim != 2."
        )
    elif (xyz.shape[0] == 3) and (xyz.shape[1] == 3):
        raise ValueError(
            "HI_disk_size: cannot calculate disk size for input with" " shape (3, 3)."
        )
    elif xyz.shape[1] == 3:
        xyz = xyz.T
        transposed = True
    elif (xyz.shape[0] != 3) and (xyz.shape[1] != 3):
        raise ValueError(
            "HI_disk_size: coordinate array shape "
            + str(xyz.shape)
            + " invalid (one dim must be 3)."
        )
    
    r = np.sqrt(np.sum(np.power(xyz, 2), axis=0))
    rsort = np.argsort(r)
    rsort = np.argsort(r, kind="quicksort")
    r = r[rsort].to('kpc')
    #print(r)
    #HI_mask = (m_HI / m_HI.max()).value > 0.3
    m_HI = m_HI[rsort]#[HI_mask]
    mcumul = np.cumsum(m_HI) / np.sum(m_HI)
    
    for j in range(len(mcumul)):
        if mcumul[j] >= frac:
            R_last = r[j] #.to('kpc')
            break
   
    #fig = plt.figure()
    #plt.plot(r, m_HI)
    #plt.scatter(R_last, 0)
    #plt.ylabel('m HI')
    #plt.xlabel('r')
    #plt.xlim([0,R_last*1.3])
    #fig.savefig(path_plots+'HI_structure_'+str(i)+'.png')
    #plt.close()
    return R_last


for simulation_type in simulation_types[6:]:
    print(simulation_type) 
    snapshotfile, halo_catalogue, catalogue_files, path_results, path_plots = manage_paths(simulation_type=simulation_type)
    
    sample_indices=np.load(path_results+'sample_indices.npy')
    
    martini_axes_Npix = []

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
    
        xyz = sg.gas.coordinates
        m_HI = sg.gas.atomic_hydrogen_masses

        R_last = HI_disk_size(xyz, m_HI, frac=0.9, save=None)
        #print(i, R_last)
        angular_size = (R_last*0.001/4)*(180/np.pi)*3600 #in arcsec
        martini_axes_Npix = []
        if angular_size/10 < 128:
            Npix = 128
        elif (angular_size/10 > 256) or (i in sample_indices[:4]):
            Npix = 512
        else:
            Npix = 256
        martini_axes_Npix.append(Npix)
        print('Npix: ', Npix)
        np.save(path_results+'R_last/R_last'+str(i)+'.npy', R_last)
        
#np.save(path_results+'martini_axes_Npix.npy', martini_axes_Npix)
print(martini_axes_Npix)
index_large_Npix = np.argwhere(martini_axes_Npix > 128)
print(index_large_Npix)