import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob
import scipy.integrate as integrate

plt.style.use("paper")

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(5.5, 5.5))

d = nc.Dataset("res_512/data.nc", "r")

pc = ax1.pcolormesh(d["x"], d["z"], d["TS"][-1,:,0,:])
pc.set_rasterized(True)
cb = plt.colorbar(pc, ax=ax1)

cb.set_label("$T$")

ax1.set_xlabel("$x$")
ax1.set_ylabel("$z$")

u = d["U"][-1,:,0,:]
u_spec = np.fft.fft2(u)
w = d["W"][-1,:,0,:]
print(u_spec.shape, w.shape)
w_spec = np.fft.fft2(w)
kx = np.fft.fftfreq(3 * 512) * 2.0 * np.pi / 10.0
kz = np.fft.fftfreq(3 * 512) * 2.0 * np.pi / 10.0

kz = 1j * kz[:,np.newaxis] - (200.0 * d["time"][-1] + d["x_remaps"][-1]) * kx[np.newaxis,:] 
kx = 0.0 * kz[:,:] + 1j * kx[np.newaxis,:]

divu_spec = 0.0 * kx
divu_spec[:,:] = u_spec * kx + w_spec * kz
print(divu_spec.shape)
divu = np.fft.ifft2(divu_spec) / 1024 / 1024

print(np.max(np.real(divu)))
pc = ax2.pcolormesh(d["x"], d["z"], np.real(divu))
pc.set_rasterized(True)
cb = plt.colorbar(pc, ax=ax2)

cb.set_label(r"$\nabla\cdot\mathbf{u}$")

ax1.set_xlabel("$x$")
ax1.set_ylabel("$z$")

labels = {"res_64": 64, "res_128": 128, "res_256": 256, "res_512": 512, "res_1024": 1024}

for dir in ["res_64", "res_128", "res_256", "res_512"]:
    d = nc.Dataset(dir + "/data.nc", "r")

    da = (d["z"][1] - d["z"][0]) * (d["x"][1] - d["x"][0])

    d = np.genfromtxt(dir + "/diag.csv", delimiter=",", names=True)
    en = 0.5 * (d["rms_u"]**2 + d["rms_v"]**2 + d["rms_w"]**2 + d["rms_t"]**2)
    prod = 200.0 * d["prod_uw"]

    ax3.plot(d["time"][1:], en[1:])
    if dir == "res_512":
        ax3.plot(d["time"][1:], en[1] - integrate.cumtrapz(prod, d["time"]), ls="--", c="k", label="$E_0-\int P_{\mathrm{sh}}dt$")
    ax4.plot(d["time"][1:], np.abs(en[1:] - en[1] + integrate.cumtrapz(prod, d["time"])), label="$N_x=%i$" % labels[dir])
    # ax3.plot(d["time"][1:], np.abs(en[1:] - en[1] + integrate.cumtrapz(prod, d["time"])), label=dir)
    # plt.plot(d["time"][1:], np.abs(en[1:] - en_interp(d["time"][1:])), label=dir)

ax3.set_xlabel("$t$")
ax3.set_ylabel("$E$")
ax3.legend(loc="lower right")
ax4.set_xlabel("$t$")
ax4.set_ylabel("$E - E_0 + \int P_{\mathrm{sh}}dt$")
ax4.set_yscale("log")
ax4.legend()
ax4.set_ylim(1e-2, 650)
plt.tight_layout()

plt.savefig("kh_test.pdf", dpi=300)
plt.show()
