import glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import interpolate

d = nc.Dataset("res_512/data.nc", "r")

da = (d["z"][1] - d["z"][0]) * (d["x"][1] - d["x"][0])

# en = 0.0 * d["time"][:]
# for i in range(len(d["time"])):
#     u = d["US"][i][:,0,:]
#     w = d["WS"][i][:,0,:]
#     T = d["TS"][i][:,0,:]
#     en[i] = np.sum(0.5*(u**2 + w**2 + 10.0*T**2)) / np.sum(u==u)

# plt.plot(d["time"][1:], np.abs(en[1:]))

d = np.genfromtxt("res_512/diag.csv", delimiter=",", names=True)

en = 0.5 * (d["rms_u"]**2 + d["rms_v"]**2 + d["rms_w"]**2 + 10.0 * d["rms_t"]**2)
en_interp = interpolate.interp1d(d["time"], en, "cubic", fill_value="extrapolate")
prod = 200.0 * d["prod_uw"]

# plt.plot(d["time"][1:], (np.abs(en[1:] - en[1] + integrate.cumtrapz(prod, d["time"]))))
# plt.plot(d["time"][1:], en[1:])

for dir in glob.glob("res_*"):
    d = nc.Dataset(dir + "/data.nc", "r")

    da = (d["z"][1] - d["z"][0]) * (d["x"][1] - d["x"][0])

    d = np.genfromtxt(dir + "/diag.csv", delimiter=",", names=True)
    en = 0.5 * (d["rms_u"]**2 + d["rms_v"]**2 + d["rms_w"]**2 + 10.0 * d["rms_t"]**2)
    prod = 200.0 * d["prod_uw"]

    plt.plot(d["time"][1:], np.abs(en[1:] - en[1]), label=dir)
    plt.plot(d["time"][1:], np.abs(integrate.cumtrapz(prod, d["time"])), label=dir)
    plt.plot(d["time"][1:], np.abs(en[1:] - en[1] + integrate.cumtrapz(prod, d["time"])), label=dir)
    # plt.plot(d["time"][1:], np.abs(en[1:] - en_interp(d["time"][1:])), label=dir)

plt.yscale("log")
plt.legend()
plt.show()

# plt.pcolor(d["x"], d["z"], T)
# plt.show()

# print(d["rms_vel"]**2 + 10.0 * d["rms_t"]**2 + )
