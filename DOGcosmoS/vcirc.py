#import DOGcosmoS
import numpy as np
import unyt as U
#from unyt import G
import matplotlib as plt
from velociraptor import load as load_catalogue
from swiftgalaxy import SWIFTGalaxy, Velociraptor
from sample_selection import read_paths_from_config

def vcirc_particle_type(xyz, m):
    """Calculates the circular velocity contributed by the particle type of interest (e.g. stars, gas or dark matter).
    
    Arguments:
        xyz :  Input array of dimension 3xN containing particle positions.
        m :    Input array of dimension 1xN containing particle masses.
             
    Returns:

        vcirc, r, m   
    """
    xyz = xyz.T
    r = np.sqrt(np.sum(np.power(xyz, 2), axis=0)).to('kpc')
    rsort = np.argsort(r, kind="quicksort")
    r = r[rsort]
    m = m[rsort]#.to('kg')
    m_encl = np.cumsum(m)
    vcirc = np.sqrt(U.G * m_encl / r).to('km/s')   
    return vcirc, r, m


def vcirc_total(halo_index, save=None, plot=None):
    """Calculates the total circular velocity.
    
    Arguments:
        halo_index: Halo index from VELOCIRaptor. Should be a single index, not an array.
        save:       If True, circular velocity profile (r,vcirc) will be saved to your output path as numpy array. Default is False
        plot:       If True, a plot showing the contribution of the circular velocity from the different particle types (i.e. stars, gas and dark matter) will be created and saved at your output path.

    Returns:
        r_all, vcirc_all
    """
    snapshotfile, halo_catalogue_filebase, halo_catalogue_properties, output_path = read_paths_from_config(config_file='specify_your_paths.ini')
    
    print('halo index', halo_index)
    
    sg = SWIFTGalaxy(
           snapshotfile,    
             Velociraptor(
             velociraptor_files=halo_catalogue_properties,
             halo_index=halo_index,
             centre_type='minpot',
             extra_mask='bound_only',
             ),
        )

    
    xyz_stars = sg.stars.coordinates
    xyz_gas = sg.gas.coordinates
    xyz_dm = sg.dark_matter.coordinates
    
    print(xyz_stars)
    print(xyz_gas)
    
    m_stars = sg.stars.masses
    m_gas = sg.gas.masses
    m_dm = sg.dark_matter.masses
    
    vcirc_stars, r_stars, m_stars = vcirc_particle_type(xyz_stars, m_stars)
    vcirc_gas, r_gas, m_gas = vcirc_particle_type(xyz_gas, m_gas)
    vcirc_dm, r_dm, m_dm = vcirc_particle_type(xyz_dm, m_dm)

    r_all = np.concatenate([r_stars, r_gas, r_dm[1:]], axis=0) # exclude innermost part because of peak behaviour
    m_all = np.concatenate([m_stars, m_gas, m_dm[1:]], axis=0)
    rsort_all = np.argsort(r_all, kind="quicksort")
    r_all = r_all[rsort_all] *U.kpc                # reintroduce lost units
    m_all = m_all[rsort_all] *10000000000*U.Msun
    m_encl_all = np.cumsum(m_all)
    vcirc_all = np.sqrt(U.G * m_encl_all / r_all).to('km/s')
    
    if save is not None:
        vcirc_profile = np.concatenate([r_all[:,np.newaxis], vcirc_all[:,np.newaxis]], axis=1)
        np.save(output_path+'r_vcirc_profile'+str(halo_index)+'.npy', vcirc_profile)
    
    if plot is not None:   
        R = sg.halo_finder.radii.r_halfmass.to('kpc').value

        fig = plt.figure()
        plt.plot(r_stars, vcirc_stars, color='tab:orange')
        plt.plot(r_gas, vcirc_gas, color='tab:cyan')
        plt.plot(r_dm[1:], vcirc_dm[1:], color='tab:purple') # exclude innermost part because of peak behaviour
        plt.plot(r_all, vcirc_all, 'k')
        plt.xlim([0,R])
        plt.title('Circular velocity')
        plt.xlabel('r [kpc]')
        plt.ylabel('v_circ [km/s]')
        plt.legend(['stars', 'gas', 'dark matter', 'total'])

        fig.savefig(output_path+'vcirc_'+str(halo_index)+'.png')
        plt.close()
        
    return r_all, vcirc_all
    
       



