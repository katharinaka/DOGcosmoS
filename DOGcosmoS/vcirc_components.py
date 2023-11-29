import numpy as np
from swiftgalaxy import Velociraptor, SWIFTGalaxy
import unyt as U
from velociraptor import load as load_catalogue
import matplotlib.pyplot as plt
from manage_paths import manage_paths
#from vrot import vrot

simulation_types = ["COLIBRE_v2022_02", "COLIBRE_v2023_03", "COLIBRE_v2023_08", "COLIBRE_v2023_08_lowres", "single-phase", "WDM", "ZOOM_Jemima"]

def vcirc_particle_type(xyz, m):    
    xyz = xyz.T
    r = np.sqrt(np.sum(np.power(xyz, 2), axis=0)).to('kpc')
    rsort = np.argsort(r, kind="quicksort")
    r = r[rsort]
    m = m[rsort]#.to('kg')
    m_encl = np.cumsum(m)
    vcirc = np.sqrt(U.G * m_encl / r).to('km/s')   
    return vcirc, r, m

for simulation_type in simulation_types:
    
    snapshotfile, halo_catalogue, catalogue_files, path_results, path_plots = manage_paths(simulation_type=simulation_type)

    sample_indices = np.load(path_results+'sample_indices.npy')
    print(simulation_type)

    for index in sample_indices:

        sg = SWIFTGalaxy(
           snapshotfile,    
             Velociraptor(
             velociraptor_files=catalogue_files,
             halo_index=index,
             centre_type='minpot',
             extra_mask='bound_only',
             ),
        )

        xyz_stars = sg.stars.coordinates
        xyz_gas = sg.gas.coordinates
        xyz_dm = sg.dark_matter.coordinates
    
        m_stars = sg.stars.masses
        m_gas = sg.gas.masses
        m_dm = sg.dark_matter.masses

        vcirc_stars, r_stars, m_stars = vcirc_particle_type(xyz_stars, m_stars)
        vcirc_gas, r_gas, m_gas = vcirc_particle_type(xyz_gas, m_gas)
        vcirc_dm, r_dm, m_dm = vcirc_particle_type(xyz_dm, m_dm)
        
        r_baryonic = np.concatenate([r_stars, r_gas], axis=0)
        m_baryonic = np.concatenate([m_stars, m_gas], axis=0)
        rsort_bar = np.argsort(r_baryonic, kind="quicksort")
        r_baryonic = r_baryonic[rsort_bar]*U.kpc                # reintroduce lost units
        m_baryonic = m_baryonic[rsort_bar]*10000000000*U.Msun
        m_encl_bar = np.cumsum(m_baryonic)
        vcirc_baryonic = np.sqrt(U.G * m_encl_bar / r_baryonic).to('km/s')
        
        vcirc_bar_profile = np.concatenate([r_baryonic[:,np.newaxis], vcirc_baryonic[:,np.newaxis]], axis=1)
        np.save(path_results+'vcirc_bar_profiles/vcirc_bar_r_profile'+str(index)+'.npy', vcirc_bar_profile)

        r_all = np.concatenate([r_stars, r_gas, r_dm[1:]], axis=0) # exclude innermost part because of peak behaviour
        m_all = np.concatenate([m_stars, m_gas, m_dm[1:]], axis=0)
        rsort_all = np.argsort(r_all, kind="quicksort")
        r_all = r_all[rsort_all]*U.kpc                # reintroduce lost units
        m_all = m_all[rsort_all]*10000000000*U.Msun
        m_encl_all = np.cumsum(m_all)
        vcirc_all = np.sqrt(U.G * m_encl_all / r_all).to('km/s')
    
        #vcirc_profile = np.concatenate([r_all[:,np.newaxis], vcirc_all[:,np.newaxis]], axis=1)
        #np.save(path_results+'vcirc_profiles/vcirc_r_profile'+str(index)+'.npy', vcirc_profile)


        R = sg.halo_finder.radii.r_halfmass.to('kpc').value

       # fig = plt.figure()
       # plt.plot(r_stars, vcirc_stars, color='tab:orange')
       # plt.plot(r_gas, vcirc_gas, color='tab:cyan')
       # plt.plot(r_dm[1:], vcirc_dm[1:], color='tab:purple') # exclude innermost part because of peak behaviour
       # plt.plot(r_all, vcirc_all, 'k')
       # plt.xlim([0,R])
       # plt.title('Circular velocity')
       # plt.xlabel('r [kpc]')
       # plt.ylabel('v_circ [km/s]')
       # plt.legend(['stars', 'gas', 'dark matter', 'total'])

       # fig.savefig(path_plots+'vcirc_'+str(index)+'.png')
