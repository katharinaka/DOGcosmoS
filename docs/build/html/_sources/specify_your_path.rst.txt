.. _specify_your_paths:

==============================
Specify Your Paths Configuration
==============================

Before diving into DOGcosmoS, you need to specify the paths to the simulation snapshot you wish to analyse, the according halo catalogue and your preferred output directory in the configuration file 'specify_your_paths.ini'. The file is structured as follows: 

.. code-block:: ini

    [INPUT]
    snapshotfile_path = /path/to/snapshotfile
    halo_catalogue_path = /path/to/halo_catalogue

    [OUTPUT]
    output_directory_path = /path/to/output_directory

Example
-------------
When working with the provided snapshot, this would be your input:


