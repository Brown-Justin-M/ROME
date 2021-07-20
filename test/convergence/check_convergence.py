import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

d = nc.Dataset("simdat01.nc", "r")

l1norm = d["time"][:] * 0.0

t0 = (100.0 / 30.0)**2 / 4.0
z = d["z"][:]
# analytic = np.sqrt(t0 / (t0 + d["time"][0])) * np.exp(-z**2 / (4 * (t0 + d["time"][0])))
analytic = np.sin(8.0 * np.pi * z / d["Gammaz"][:])
# a = np.mean((d["Temp"][0,:,0,0] - analytic))
a = 0.0
print(np.max(d["z"]))
for i in range(len(d["time"])):
    # line = plt.plot(d["z"], d["Temp"][i,:,0,0] - a)[0]
    # analytic = np.sqrt(t0 / (t0 + d["time"][i])) * np.exp(-z**2 / (4 * (t0 + d["time"][i])))
    analytic = np.exp(-(8.0 * np.pi / d["Gammaz"][:])**2 * d["time"][i]) * np.sin(8.0 * np.pi * (z - 1.0*d["time"][i]) / d["Gammaz"][:])
    # plt.plot(d["z"], analytic, c=line.get_color(), ls="--")
    # plt.plot(d["z"], d["Temp"][i,:,0,0] - analytic,ls="--")

    l1norm[i] = np.sqrt(np.mean(((d["Temp"][i,:,0,0] - a) - analytic)**2) / np.mean(analytic**2))

plt.plot(d["time"], l1norm)
print("%15.12f" % np.log10(l1norm[-1]))
plt.show()