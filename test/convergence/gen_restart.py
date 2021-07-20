import f90nml
import netCDF4 as nc
import numpy as np

file = open("parameters", "r")
nl = f90nml.read(file)

grid = nl["grid"]
params = nl["parameters"]["params"]
io = nl["io"]["controls"]

nx = grid["x_modes"]
ny = 0
nz = grid["z_modes"]

d = nc.Dataset("infile.nc", "w", format="NETCDF3_64BIT")

xid = d.createDimension("kx", nx + 1)
yid = d.createDimension("ky", max(2 * ny, 1))
zid = d.createDimension("kz", 2 * nz)

dims = ("ky", "kx", "kz")

d.createVariable("kx", "f8", "kx")
d.createVariable("ky", "f8", "ky")
d.createVariable("kz", "f8", "kz")

vars = ["T", "C", "U", "V", "W", "TI", "CI", "UI", "VI", "WI"]
for var in vars:
    d.createVariable(var, "f8", dims)
    d[var][:] = 0.0

gammaz = grid["z_length"]

kz = np.fft.fftfreq(2 * nz) * 2 * nz * 2.0 * np.pi / gammaz

d["T"][0,0,:] = np.exp(-(kz[np.newaxis,np.newaxis,:])**2 / (2.0e2)**2)
# d["T"][:] = np.sin(8.0 * np.pi * z / gammaz)

d.close()
