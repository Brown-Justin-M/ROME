import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("nz", type=int)
parser.add_argument("dt", type=float)

args = parser.parse_args()

nz = args.nz
dt = args.dt

steps = int(0.01 / dt)

params = "&grid\nx_modes = 2\nz_modes = %i\nz_length = 1.0\nz_vel = 100.0\n/\n\n" % (nz)

params = params + "&parameters\nparams%%max_time = 0.01\nparams%%start_dt = %e\nparams%%max_dt = %e\n" % (dt,dt)

params = params + """params%diff_vel = 0.0
params%strat_t = 0.0
params%buoy_t = 0.0
params%diff_t = 1.0
params%strat_c = 0.0
params%buoy_c = 0.0
params%diff_c = 0.0
/

&io
"""

params = params + "controls%%steps_data = %i\n" % (int(steps / 10))

params = params + """controls%steps_restart = 100
controls%steps_diag = 1
controls%steps_profile = 100
controls%file_input = "infile.nc"
/
"""

file = open("parameters", "w")

file.write(params)

file.close()
