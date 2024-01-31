# DOGcosmoS

Dynamics Of Galaxies in Cosmological Simulations

# Installation
Having cloned this project from github, you can move into the root directory and install via
pip -m install --editable .

Please make sure that you have the requirements installed before using this package:

swiftgalaxy == 1.0.0
velocyraptor
numpy == 1.23.5
unyt == 2.9.3

Since the rewuirements require specific versions to work well within this package, it is recommended to use a virtual environment.

# Usage

For anyone  who would like to analyse the dynamics of galaxies in a cosmological simulation box.

A quick-start guide:
First you need to specify the paths to the snapshotfile and halo catalogue in the file 'specify_your_paths.ini'.

In your script, you need to import sample_selection and other functions if needed. For example, if you wish to plot the rotation curve diversity for a galaxy sample with vmax between 60 and 120 km/s your script could look like this:

from sample_selection import read_paths_from_config, select_sample
from vcirc imort vcirc_total
from rotation_curve_shapes import quantify_rotation_curve_shapes
from plot_rotation_curve_diversity import plot_vfid_vmax

plot_vfid_vmax(60,120)


# Contribution
This package is open source, you can contribute via github.

# License
This project is licensed under the GNU general licence.

