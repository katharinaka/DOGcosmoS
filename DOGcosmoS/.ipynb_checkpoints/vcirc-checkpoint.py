import numpy as np
import unyt as U
from velociraptor import load as load_catalogue

def vcirc_particle_type(xyz, m, save=False):
    """Calculates the circular velocity contributed by the particle type of interest (e.g. stars, gas or dark matter).
    
    Arguments:
             xyz:  Input array of dimension 3xN containing particle positions.
             m:    Input array of dimension 1xN containing particle masses.
             save: If True, circular velocity profile (r,vcirc) will be saved to your output path as numpy array. Default is False.
    """
    xyz = xyz.T
    r = np.sqrt(np.sum(np.power(xyz, 2), axis=0)).to('kpc')
    rsort = np.argsort(r, kind="quicksort")
    r = r[rsort]
    m = m[rsort]#.to('kg')
    m_encl = np.cumsum(m)
    vcirc = np.sqrt(U.G * m_encl / r).to('km/s')   
    return vcirc, r, m

def vcirc_total(xyz_stars, xyz_gas, xyz_dm,
                m_stars, m_gas, m_dm,
                save=False):
    """Calculates the total circular velocity.
    
    Arguments:
             xyz_stars:  Input array of dimension 3xN containing the positions of star particles.
             xyz_gas:    Definition analogous to xyz_stars.
             xyz_dm:     Definition analogous to xyz_stars.
             m_star:     Input array of dimension 1xN containing the masses of star particles.
             m_gas:      Definition analogous to m_star.
             m_dm:       Definition analogous to m_star.
             save:       If True, circular velocity profile (r,vcirc) will be saved to your output path as numpy array. Default is False
    """
    vcirc_stars, r_stars, m_stars = vcirc_particle_type(xyz_stars, m_stars, save=False)
    vcirc_gas, r_gas, m_gas = vcirc_particle_type(xyz_gas, m_gas, save=False)
    vcirc_dm, r_dm, m_dm = vcirc_particle_type(xyz_dm, m_dm, save=False)
        
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
    
    return r_all, vcirc_all
    
        #vcirc_profile = np.concatenate([r_all[:,np.newaxis], vcirc_all[:,np.newaxis]], axis=1)
        #np.save(path_results+'vcirc_profiles/vcirc_r_profile'+str(index)+'.npy', vcirc_profile)


    #R = sg.halo_finder.radii.r_halfmass.to('kpc').value

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
