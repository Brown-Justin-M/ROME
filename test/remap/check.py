import numpy as np
import netCDF4 as nc
import f90nml

import matplotlib.pyplot as plt

file = open("parameters", "r")
nl = f90nml.read(file)

grid = nl["grid"]

d = nc.Dataset("data.nc","r")

vals = []
for i in [0,len(d["time"]) - 1]:
    t = d["TS"][i][:,0,:]
    z = d["z"][:]
    z = z[:,np.newaxis] + 0.0 * t
    x = d["x"][:]
    x = (x[np.newaxis,:] - grid["x_shear"] * d["time"][i] * z) % 100.0
    z = z - grid["z_length"] / 2

    # r = np.sqrt((x - 50.0)**2 + (z - 50.0)**2)
    # analytic = np.exp(-r**2 / 100.0)

    fac = np.cos(np.pi * z / 50.0)**2
    mask = np.logical_and(z >= -25.0, z <= 25.0)
    fac[np.logical_not(mask)] = 0.0
    analytic =  np.tanh((z + 12.5 * np.sin(2.0 * np.pi * x / 100.0)) / 1.0) * fac

    vals.append(np.sum(np.abs(t[mask] - analytic[mask])) / np.sum(mask))

    # pc = plt.pcolormesh(d["x"], d["z"], analytic)
    # cb = plt.colorbar(pc)
    # plt.show()

print(grid["x_modes"], grid["z_modes"], grid["z_length"] / grid["x_length"], vals[0], vals[1])