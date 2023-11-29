import numpy as np
from swiftgalaxy import Velociraptor, SWIFTGalaxy
import unyt as U
from velociraptor import load as load_catalogue
import matplotlib.pyplot as plt
from scipy.stats import tmean, binned_statistic

snapshotfile = "/gpfs/data/fs71897/kkain/snapshot_data/colibre_0024.hdf5"

halo_catalogue = "/gpfs/data/fs71897/kkain/snapshot_data/halos_0024"
halo_catalogue_file = "/gpfs/data/fs71897/kkain/snapshot_data/halos_0024.properties"

path_plots = '/gpfs/data/fs71897/kkain/Plots/'

from rotmat import L_align



def azimuthal_velocity(xyz, vxyz, m, frac=0.5, saverot=None):
    
    xyz = sg.gas.coordinates
    vxyz = sg.gas.velocities
    
    Lx = sg.halo_finder.angular_momentum.lx_200m_gas
    Ly = sg.halo_finder.angular_momentum.ly_200m_gas
    Lz = sg.halo_finder.angular_momentum.lz_200m_gas

    