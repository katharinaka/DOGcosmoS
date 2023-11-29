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

def NFW_vcirc(vvir, r, rvir, c):
    x = r/rvir
    fc = np.log(1+c) - (c / (1+c)) # or +?
    fcx = np.log(1+ c*x) - (c*x / (1+ c*x))
    vcirc = vvir * np.sqrt(fcx / (x*fc))
    return vcirc

def estimate_Delta_vir(Omega_m):
    y = Omega_m - 1
    Delta_vir = (18 * np.pi**2 +82*y - 39*y**2)/(y+1)
    print(Delta_vir*Omega_m)
    return Delta_vir

def estimate_rvir(Mvir, Delta_vir, Omega_m, h, z):
    rvir = 163 *h**-1  *U.kpc * (Mvir/ (10**12 *h**-1 *U.Msun) )**(1/3) * (Delta_vir*Omega_m/200)**-(1/3)* (1+z)**-1 #( 1/Omega_m = Delta_vir/200  (?))
    return rvir

def estimate_vvir(Mvir, Delta_vir, Omega_m, h, z):
    vvir = 163  *(U.km / U.s) * (Mvir/ (10**12 * h**-1 *U.Msun) )**(1/3) * (Delta_vir*Omega_m/200)**-(1/6)* np.sqrt(1+z)
    return vvir

def estimate_vmax(vvir, c):
    fc = np.log(1+c) - (c / (1+c))
    vmax = 0.465 * vvir * np.sqrt(c/fc)
    return vmax


#def NFW_concentration(M, z):
    #z = 0
#    alpha = 1.62774 - 0.2458 * (1. + z) + 0.01716 * (1. + z)**2
#    beta = 1.66079 + 0.00359 * (1. + z) - 1.6901 *(1. + z)**0.00417
#    gamma = -0.02049 + 0.0253 * (1. + z)**(-0.1044)#

#    log_c = alpha + beta * np.log10(M/U.Msun) * (1. + gamma * np.log10(M/U.Msun)**2 )
#    c = np.power(10, log_c)
#    return c

def NFW_concentration(M,z):
    alpha = 1.7543 - 0.2766 * (1 + z) + 0.02039 * (1 + z)**2
    beta = 0.2753 + 0.00351 * (1 + z) - 0.3038 * (1 + z)**0.0269
    gamma = -0.01537 + 0.02102 * (1 + z)**(-0.1475)
    log_c = alpha + beta * np.log10(M/U.Msun) * (1 + gamma * (np.log10(M/U.Msun)**2))
    c = 10**log_c
    return c

mass_range = np.linspace(1e5, 1e15, 100)*U.Msun.to('Msun')
print('log mass range :', np.log10(mass_range/U.Msun))
#Planck1:
#if simulation_type == "COLIBRE_v2022_02" or simulation_type == "single-phase":
#print(repeat)
#    Omega_m = 0.307
#    h = 0.6777
#else:    
Omega_m = 0.30
h = 0.681

c = NFW_concentration(M=mass_range, z=0)

sigma_logc = 0.1
c_upper = c + c * 0.259 # (10**0.1-1)
print('c_upper: ', c_upper)
c_lower = c - c * 0.2589 #(10**0.1-1)
print('c_lower: ', c_lower)
#c = np.load('NFW_concen
#import commah
#output = commah.run('WMAP5',zi=0.,Mi=mass_range ,z=0.)
#c = output['c'].flatten()

#c = np.load('NFW_concentration_commah.npy')
print('c : ', c)
z=0

#WMAP9:
#Omega_m = 0.282
#h = 0.70
#WMAP1:
#Omega_m = 0.25
#h = 0.73
#WMAP5:
#Omega_m = 0.258
#h = 0.72
Delta_vir = estimate_Delta_vir(Omega_m)
rvir_NFW = estimate_rvir(mass_range, Delta_vir, Omega_m, h, z)
vvir_NFW = estimate_vvir(mass_range, Delta_vir, Omega_m, h, z)
vmax_NFW = estimate_vmax(vvir_NFW, c)
print('vmax: ', vmax_NFW)
print('rvir: ', rvir_NFW)
rfid_NFW = 2 * vmax_NFW/ (70 *U.km / U.s) * U.kpc
print('rfid: ', rfid_NFW)
vfid_NFW = NFW_vcirc(vvir_NFW, rfid_NFW, rvir_NFW, c)
print('vfid: ', vfid_NFW)
np.save('vmax_vfid_NFW.npy', np.concatenate([vmax_NFW[:,np.newaxis], vfid_NFW[:,np.newaxis]], axis=1))

#vmax_upper = estimate_vmax(vvir_NFW, c_upper)
#print('vmax upper: ', vmax_upper)
#print('vmax up with clow: ',estimate_vmax(vvir_NFW, c_lower))
#rfid_upper = 2 * vmax_upper / (70 * U.km / U.s) *U.kpc
#vfid_upper = NFW_vcirc(vvir_NFW, rfid_upper, rvir_NFW, c_upper)
#vmax_lower = estimate_vmax(vvir_NFW, c_lower)    
#rfid_lower = 2 * vmax_lower / (70 * U.km / U.s) *U.kpc
#vfid_lower = NFW_vcirc(vvir_NFW, rfid_lower, rvir_NFW, c_lower)

simulation_types = ["COLIBRE_v2022_02", "COLIBRE_v2023_03", "COLIBRE_v2023_08", "COLIBRE_v2023_08_lowres", "single-phase", "WDM", "ZOOM_Jemima"]

#comparefig = plt.figure()
#plot_line = np.linspace(59,121)
#plt.plot(plot_line, plot_line, "k:")
#plt.plot(vmax_NFW, vfid_NFW, color="grey")

for simulation_type in simulation_types[6:]:
    print(simulation_type) 
    snapshotfile, halo_catalogue, catalogue_files, path_results, path_plots = manage_paths(simulation_type=simulation_type)
    
    sample_indices=np.load(path_results+'sample_indices.npy')
    
    vmax_sample = np.zeros(len(sample_indices))
    r_fid_sample = np.zeros(len(sample_indices))
    v_fid_sample = np.zeros(len(sample_indices))
    eta_rot_sample = np.zeros(len(sample_indices))

    for i in range(len(sample_indices)):
        index = sample_indices[i]
    
        vcirc_profile = np.load(path_results+'vcirc_profiles/vcirc_r_profile'+str(index)+'.npy')
        r = vcirc_profile[:,0]
        vcirc = vcirc_profile[:,1]
    
        vmax = np.max(vcirc)
        r_fid = 2 * vmax/70   #/70km/s kpc

        # Perform spline interpolation
        spl = splrep(r, vcirc)

        # Evaluate the spline at desired points
        v_fid = splev(r_fid, spl)
    
        vmax_sample[i] = vmax
        r_fid_sample[i] = r_fid
        v_fid_sample[i] = v_fid
        eta_rot_sample[i] = v_fid/vmax

    vmax_vfid = np.concatenate([vmax_sample[:,np.newaxis], v_fid_sample[:,np.newaxis]], axis=1)
    np.save(path_results+'vmax_vfid.npy', vmax_vfid)

    coeff_1d = np.polyfit(vmax_sample, v_fid_sample, 1)
    coeff_2d = np.polyfit(vmax_sample, v_fid_sample, 2)
    fit_1d = np.poly1d(coeff_1d)
    fit_2d = np.poly1d(coeff_2d)


#fig2 = plt.figure()
#plt.plot(np.log10(mass_range/U.Msun), c)
# #plt.xscale('log')
#plt.yscale('log')
#fig2.savefig(path_plots+'M-c_relation.png')
#plt.close()

    vfid_upper = vfid_NFW * 10**0.1 # (10**0.1-1)
    print('vfid_upper: ', vfid_upper)
    vfid_lower = vfid_NFW * 10**-0.1 #(10**0.1-1)
    print('vfid_lower: ', vfid_lower)

    #flag_bad_indices = [26, 30, 34, 35, 36, 37, 40, 42, 43, 48, 50, 51, 54, 55, 56, 57, 66]
    
    fig = plt.figure()
    plot_line = np.linspace(20,150)
    plt.plot(plot_line, plot_line, 'k:')
    #plt.plot(bins[:-1], running_medians, label='Adaptive Running Median', color='red')
    plt.plot(vmax_NFW, vfid_NFW, color="grey")
    plt.fill_between(vmax_NFW, vfid_lower, vfid_upper, color='gray', alpha=0.3, label='_nolegend_')
    data = plt.scatter(vmax_sample, v_fid_sample, c=eta_rot_sample, vmin=0, vmax=1, cmap='viridis')
    #for i in range(len(sample_indices)):
       #index = sample_indices[i]
       # if index in flag_bad_indices:
       #     color = "red"
       # else:
       #     color = "black"
       # plt.scatter(vmax_sample[i], v_fid_sample[i], alpha=0)
       # plt.annotate(index, (vmax_sample[i], v_fid_sample[i]), color=color)
    
    plt.xlabel('$V_{max}$', fontsize=18)
    plt.ylabel('$V_{fid}$', fontsize=18)
    #plt.xlim([59, 121])
    #plt.ylim([25, 120])
    plt.xlim([20, 150])
    plt.ylim([3, 120])
    
    if simulation_type == "single-phase":
        legend = 'single-phase ISM, $V_{circ}$'
    else:
        legend = 'multi-phase ISM, $V_{circ}$'
    plt.legend(['1:1', 'NFW', legend]) 
    plt.xscale('log')
    plt.yscale('log')
    #x_ticks = [60, 70, 80, 90, 100, 110, 120]  # These are the true numbers
    #x_tick_labels = ['60', '70', '80', '90', '100', '110', '120']  # Custom labels
    #plt.xticks(x_ticks, x_tick_labels, fontsize=16)
    #y_ticks = [30, 40, 50, 60, 70, 80, 90, 100]  # These are the true numbers
    #y_tick_labels = ['30', '40', '50', '60', '70', '80', '90', '100']  # Custom labels
    #plt.yticks(y_ticks, y_tick_labels, fontsize=16)
    
    x_ticks = [20, 30, 40, 50, 60, 80, 100, 120, 150]  # These are the true numbers
    x_tick_labels = ['20', '30', '40', '50', '60', '80', '100', '120', '150']  # Custom labels
    plt.xticks(x_ticks, x_tick_labels, fontsize=16)
    y_ticks = [10, 20, 30, 50, 70, 90, 120]  # These are the true numbers
    y_tick_labels = ['10', '20', '30', '70', '50', '90', '120']  # Custom labels
    plt.yticks(y_ticks, y_tick_labels, fontsize=16)
    
    plt.plot(np.ones(117)*60, np.linspace(3,120, 117), linestyle="--", color="grey")
    plt.plot(np.ones(117)*120, np.linspace(3,120, 117), linestyle="--", color="grey")
    #ax.axis([59, 121, 25, 120])
    #ax.loglog()
    #for axis in [ax.xaxis, ax.yaxis]:
    #    formatter = ScalarFormatter()
    #    formatter.set_scientific(False)
    #    axis.set_major_formatter(formatter)
    #plt.title('vcirc, multiphase')
    cbar = fig.colorbar(data)
    cbar.set_label('$\eta_{rot} = v_{fid} / v_{max}$')
    plt.tight_layout()
    fig.savefig(path_plots+'vfid_vs_vmax_circ_fullrange.png', dpi=400)
    plt.close()



    # Perform spline interpolation
    spl = splrep(vmax_NFW, vfid_NFW)

    # Evaluate the spline at desired points
    NFW_vfid = splev(vmax_sample, spl)
    
    dev = v_fid_sample/NFW_vfid
    
    fig2 = plt.figure(figsize=(7,6))
    plt.plot(vmax_NFW, vfid_NFW/vfid_NFW, color="grey")
    plt.fill_between(np.linspace(59,121,len(vfid_NFW)), vfid_NFW/vfid_NFW*10**-0.1, vfid_NFW/vfid_NFW*10**0.1, color='gray', alpha=0.3, label='_nolegend_')
    plt.scatter(vmax_sample, dev)
    plt.xlim([59,121])
    plt.ylim([0.4,1.6])
    plt.xlabel('$V_{max}$', fontsize=18)
    plt.ylabel('$V_{fid}\  / \ V_{fid,NFW}$', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    fig2.savefig(path_plots+'deviation_from_NFW.png', dpi=150)
    plt.close()
