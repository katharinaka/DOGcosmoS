import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from swiftgalaxy import Velociraptor, SWIFTGalaxy
import unyt as U
from velociraptor import load as load_catalogue
from scipy.interpolate import splrep, splev
from scipy.stats import binned_statistic
import astropy
from manage_paths import manage_paths

simulation_types = ["COLIBRE_v2022_02", "COLIBRE_v2023_03", "COLIBRE_v2023_08", "COLIBRE_v2023_08_lowres", "single-phase", "WDM", "ZOOM_Jemima", "TANGO_CDM", "TANGO_Ref_SigmaVelDep_60", "TANGO_WeakStellarFB_SigmaVelDep_60"]

for simulation_type in simulation_types:
    print(simulation_type) 
    snapshotfile, halo_catalogue, catalogue_files, path_results, path_plots = manage_paths(simulation_type=simulation_type)
    
    sample_indices=np.load(path_results+'sample_indices.npy')
    
    vmax_sample = np.zeros(len(sample_indices))
    r_fid_sample = np.zeros(len(sample_indices))
    v_fid_sample = np.zeros(len(sample_indices))
    eta_rot_sample = np.zeros(len(sample_indices))
    eta_bar_sample = np.zeros(len(sample_indices))

    for i in range(len(sample_indices)):
        index = sample_indices[i]
    
        vcirc_profile = np.load(path_results+'vcirc_profiles/vcirc_r_profile'+str(index)+'.npy')
        r = vcirc_profile[:,0]
        vcirc = vcirc_profile[:,1]
        
        vcirc_bar_profile = np.load(path_results+'vcirc_bar_profiles/vcirc_bar_r_profile'+str(index)+'.npy')
        r_bar = vcirc_bar_profile[:,0]
        vcirc_bar = vcirc_bar_profile[:,1]
    
        vmax = np.max(vcirc)
        r_fid = 2 * vmax/70   #/70km/s kpc

        # Perform spline interpolation
        spl = splrep(r, vcirc)
        spl_bar = splrep(r_bar, vcirc_bar)

        # Evaluate the spline at desired points
        v_fid = splev(r_fid, spl)
        v_bar_fid = splev(r_fid, spl_bar)
    
        vmax_sample[i] = vmax
        r_fid_sample[i] = r_fid
        v_fid_sample[i] = v_fid
        eta_rot_sample[i] = v_fid/vmax
        eta_bar_sample[i] = (v_bar_fid/v_fid)**2
    

    np.save(path_results+'eta_rot', eta_rot_sample)
    np.save(path_results+'eta_bar', eta_bar_sample)


    fig = plt.figure(figsize=(7,6))
    #plt.plot(vmax_NFW, vfid_NFW/vfid_NFW, color="grey")
    #plt.fill_between(np.linspace(59,121,len(vfid_NFW)), vfid_NFW/vfid_NFW*10**-0.1, vfid_NFW/vfid_NFW*10**0.1, color='gray', alpha=0.3, label='_nolegend_')
    data = plt.scatter(eta_bar_sample, eta_rot_sample, c=vmax_sample, cmap='viridis', vmin=60, vmax=120)
    cbar = fig.colorbar(data)
    cbar.set_label('$V_{max}$')
    #plt.xlim([0,1])
    #plt.ylim([0.2,1])
    plt.xlim([0.08,1])
    plt.ylim([0.2,1])
    plt.xscale('log')
    plt.yscale('log')
    #plt.yticks(color='w')
    x_ticks = [0.1, 0.3, 0.7, 1.0]  # These are the true numbers
    x_tick_labels = ['0.1', '0.3', '0.7', '1.0']  # Custom labels
    plt.xticks(x_ticks, x_tick_labels, fontsize=16)
    y_ticks = [0.2, 0.3, 0.4, 0.6, 1.0]  # These are the true numbers
    y_tick_labels = ['0.2', '0.3', '0.4', '0.6', '1.0']  # Custom labels
    plt.yticks(y_ticks, y_tick_labels, fontsize=16)

    plt.xlabel('$\eta_{bar} = (V_{bar, fid}\ /\ V_{fid})^2$', fontsize=18)
    plt.ylabel('$\eta_{rot} = V_{fid}\ /\ V_{max}$', fontsize=18)
    #plt.xticks(fontsize=16)
    #plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.title(str(simulation_type))
    fig.savefig(path_plots+'etarot_vs_etabar.png')
    plt.close()
