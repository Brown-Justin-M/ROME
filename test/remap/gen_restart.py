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
zid = d.createDimension("kz", max(2 * nz, 1))

dims = ("ky", "kx", "kz")

d.createVariable("kx", "f8", "kx")
d.createVariable("ky", "f8", "ky")
d.createVariable("kz", "f8", "kz")

d.createVariable("T", "f8", dims)
d.createVariable("C", "f8", dims)
d.createVariable("U", "f8", dims)
d.createVariable("V", "f8", dims)
d.createVariable("W", "f8", dims)

d.createVariable("TI", "f8", dims)
d.createVariable("CI", "f8", dims)
d.createVariable("UI", "f8", dims)
d.createVariable("VI", "f8", dims)
d.createVariable("WI", "f8", dims)

print(nx + 1,max(2 * ny, 1),max(2 * nz, 1))

kx = np.fft.rfftfreq(2 * nx) * 2.0 * np.pi / grid["x_length"]
kz = np.fft.fftfreq(2 * nz) * 2.0 * np.pi / grid["z_length"]
kx = np.array(d["T"][:]) * 0.0 + kx[np.newaxis,:,np.newaxis]
kz = np.array(d["T"][:]) * 0.0 + kz[np.newaxis,np.newaxis,:]

kr = np.sqrt(kx**2 + kz**2)

print(np.max(kr))

# d["T"][:] = kr * np.exp(-(kr - 0.005)**2 / 1.0e-6) * np.sin(4.0 * np.arctan2(kz, kx))
d["T"][:] = np.exp(-(kr)**2 / 1.0e-7)
d["C"][:] = 0.0
d["U"][:] = 0.0
d["V"][:] = 0.0
d["W"][:] = 0.0

d["TI"][:] = 0.0
d["CI"][:] = 0.0
d["UI"][:] = 0.0
d["VI"][:] = 0.0
d["WI"][:] = 0.0

d.close()
