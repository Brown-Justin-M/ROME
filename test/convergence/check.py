import numpy as np
import netCDF4 as nc
import f90nml

import matplotlib.pyplot as plt

file = open("parameters", "r")
nl = f90nml.read(file)

grid = nl["grid"]
params = nl["parameters"]["params"]

d = nc.Dataset("data.nc","r")

gammaz = grid["z_length"]
upfactor = 2048 // grid["z_modes"]
kz = np.zeros(3 * upfactor * grid["z_modes"], dtype=np.double)
kz[:] = np.fft.fftfreq(3 * upfactor * grid["z_modes"]) * 3 * upfactor * grid["z_modes"] * 2.0 * np.pi / gammaz

init = np.exp(-(kz)**2 / (2.0e2)**2)

vals = []
for i in [0,len(d["time"]) - 1]:
# for i in range(len(d["time"])):
    t = d["T"][i][:,0,0]
    z = d["z"][:]
    zshift = d["time"][i] * grid["z_vel"]

    solution = init * np.exp(-params["diff_t"] * kz**2 * d["time"][i] + 1j * kz * zshift)
    analytic = np.fft.fft(solution)[::upfactor]

    # plt.plot(z,np.real(analytic))
    # plt.plot(z,t)
    # plt.show()

    vals.append(np.sum(np.abs(t - np.real(analytic))) / np.sum(analytic == analytic))
    # vals.append(np.abs(t[0] - analytic[0]))

print(grid["z_modes"], d["dt"][-1], vals[0], vals[1])