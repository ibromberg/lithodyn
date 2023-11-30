# Plotting from COMSOL
These are my pre- and post-processing files! 

Process for plotting velocities:
Export final topography data from comsol (in xyz); run in datagrid.py to generate a regular grid along lat/long at which to export velocities at. This generates two files - one for topography and one for 28 km below sea level. Take these files and use them in comsol's export data ("points to evaluate in: from file"). Then you can export velocity vectors (files in format of x y z u v w) at topography and 28 km bsl. These two files get uploaded to diffflow.py to create a plot of surface/28kmbsl velocities and a plot of differential velocities, with core complex data. 

topo_changes.py just plots topography for 3 different time steps. 
