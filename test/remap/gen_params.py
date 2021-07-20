import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("nx", type=int)
parser.add_argument("aspect", type=float)

args = parser.parse_args()

nx = args.nx
aspect = args.aspect

nz = int((np.ceil(nx / 2) * aspect)) * 2

x_length = 100.0
z_length = x_length * aspect

shear = 0.1 * x_length / z_length

params = "&grid\nx_modes = %i\nz_modes = %i\nx_length = 100.0\nz_length = %f\nx_shear = %f\n/\n\n" % (nx, nz, z_length, shear)

params = params + """&parameters
params%max_steps = 100
params%start_dt = 1.E-1
params%max_dt = 1.E-1
params%diff_vel = 0.0
params%strat_t = 0.0
params%buoy_t = 0.0
params%diff_t = 0.0
params%strat_c = 0.0
params%buoy_c = 0.0
params%diff_c = 0.0
/

&io
controls%steps_data = 10
controls%steps_restart = 100
controls%steps_diag = 1
controls%steps_profile = 100
controls%file_input = "infile.nc"
controls%stationary_output = .true.
/
"""

file = open("parameters", "w")

file.write(params)

file.close()
