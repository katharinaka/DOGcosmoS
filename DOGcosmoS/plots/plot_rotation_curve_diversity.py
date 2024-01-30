import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import unyt as U
from velociraptor import load as load_catalogue
from matplotlib.ticker import ScalarFormatter

#from scipy.interpolate import splrep, splev
#from scipy.stats import binned_statistic
import astropy

from sample_selection import select_sample, read_paths_from_config
from rotation_curve_shapes import quantify_rotation_curve_shapes

def NFW_vcirc(vvir, r, rvir, c):
    x = r/rvir
    fc = np.log(1+c) - (c / (1+c)) # or +?
    fcx = np.log(1+ c*x) - (c*x / (1+ c*x))
    vcirc = vvir * np.sqrt(fcx / (x*fc))
    return vcirc

def estimate_Delta_vir(Omega_m):
    y = Omega_m - 1
    Delta_vir = (18 * np.pi**2 +82*y - 39*y**2)/(y+1)
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

def plot_vfid_vmax(vmax_min_sample, vmax_max_sample, title='Rotation curve diversity', legend='simulation data', plot_NFW_ref=True, save_as='vfid_vs_vmax.png', dpi=200):
#def plot_vfid_vmax(sample_indices, title='Rotation curve diversity', legend='simulation data', plot_NFW_ref=True, save_as='vfid_vs_vmax.png', dpi=200):

    vmax_sample, vfid_sample, eta_rot_sample = quantify_rotation_curve_shapes(vmax_min_sample, vmax_max_sample, save=None)
    
    fig = plt.figure()
    plot_line = np.linspace(vmax_min_sample, vmax_max_sample)
    plt.plot(plot_line, plot_line, 'k:')
    
    if plot_NFW_ref:

        Omega_m = 0.30
        h = 0.681
        z = 0
        mass_range = np.linspace(1e5, 1e15, 100)*U.Msun.to('Msun')

        c = NFW_concentration(M=mass_range, z=0)
        Delta_vir = estimate_Delta_vir(Omega_m)
        rvir_NFW = estimate_rvir(mass_range, Delta_vir, Omega_m, h, z)
        vvir_NFW = estimate_vvir(mass_range, Delta_vir, Omega_m, h, z)
        vmax_NFW = estimate_vmax(vvir_NFW, c)
        rfid_NFW = 2 * vmax_NFW/ (70 *U.km / U.s) * U.kpc
        vfid_NFW = NFW_vcirc(vvir_NFW, rfid_NFW, rvir_NFW, c)

        vfid_upper = vfid_NFW * 10**0.1 # (10**0.1-1)
        vfid_lower = vfid_NFW * 10**-0.1 #(10**0.1-1)
    
        plt.plot(vmax_NFW, vfid_NFW, color="grey")
        plt.fill_between(vmax_NFW, vfid_lower, vfid_upper, color='gray', alpha=0.3, label='_nolegend_')
        
    data = plt.scatter(vmax_sample, vfid_sample, c=eta_rot_sample, vmin=0, vmax=1, cmap='viridis')
    
    plt.xlabel('$V_{max}$', fontsize=18)
    plt.ylabel('$V_{fid}$', fontsize=18)
    ##plt.xlim([59, 121])
    ##plt.ylim([25, 120])
    plt.xlim([20, 150])
    plt.ylim([3, 120])
    plt.title(title)
    if plot_NFW_ref:
        plt.legend(['1:1', 'NFW', legend])
    else:
        plt.legend(['1:1', legend])
    plt.xscale('log')
    plt.yscale('log')
    
    x_ticks = [20, 30, 40, 50, 60, 80, 100, 120, 150]  # These are the true numbers
    x_tick_labels = ['20', '30', '40', '50', '60', '80', '100', '120', '150']  # Custom labels
    plt.xticks(x_ticks, x_tick_labels, fontsize=16)
    y_ticks = [10, 20, 30, 50, 70, 90, 120]  # These are the true numbers
    y_tick_labels = ['10', '20', '30', '70', '50', '90', '120']  # Custom labels
    plt.yticks(y_ticks, y_tick_labels, fontsize=16)
    
    plt.plot(np.ones(117)*60, np.linspace(3,120, 117), linestyle="--", color="grey")
    findplt.plot(np.ones(117)*120, np.linspace(3,120, 117), linestyle="--", color="grey")
    cbar = fig.colorbar(data)
    cbar.set_label(r'$\eta_{rot} = v_{fid} / v_{max}$')
    plt.tight_layout()
    snapshotfile, halo_catalogue_files, output_path = read_paths_from_config()
    fig.savefig(output_path+save_as, dpi=dpi)
    plt.close()
                
    return

