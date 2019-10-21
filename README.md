# roms_sigma_coordinate
Code for data in ROMS native coordinate system (i.e. sigma coordinates)

set_depth.py calls stretching to compute the ROMS depths for the various grids used internally by ROMS.</br>
igrid=1 is the grid for the density variables (temp, salinity, etc).</br>
igrid=2 is for stream function points</br>
igrid=3 is the grid for u-velcoity points</br>
igrid=4 is the grid for v-velocity points</br>
igrid=5 is the grid for w-velocity points</br>

The compute_heat_content_MBNMS computes the heat content for a single time point (time=0) in the model output.
It is an example of how to call the various parameters and get an output.
A simple trapezoid rule is used to compute the vertical integral at a grid point.
Final output is in Giga Joules/m^2.

Ideally you would want to look at the change between times at a grid point which should be a smaller number.

