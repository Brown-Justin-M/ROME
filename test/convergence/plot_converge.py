import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

plt.style.use("paper")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.5, 3.0))

d = nc.Dataset("data.nc", "r")

gammaz = 1.0
upfactor = 2
kz = np.zeros(3 * upfactor * 1024, dtype=np.double)
kz[:] = np.fft.fftfreq(3 * upfactor * 1024) * 3 * upfactor * 1024 * 2.0 * np.pi / gammaz

init = np.exp(-(kz)**2 / (2.0e2)**2)
z = d["z"][:]

for i in range(0, len(d["time"]), 2):
    t = d["T"][i,:,0,0]
    t[2764-10:2764+10] = np.nan
    ax1.plot((d["z"][:] + 0.1) % 1.0 - 0.1, t)

    zshift = d["time"][i] * 100.0

    solution = init * np.exp(-1.0 * kz**2 * d["time"][i] + 1j * kz * zshift)
    analytic = np.fft.fft(solution)[::upfactor]
    # ax1.plot((d["z"][:] - 0.5) % 1.0 - 0.5, analytic, ls="--")

t = np.linspace(0, 0.01, 100)
zshift = 0.0 * t
maxes = 0.0 * t
for i in range(len(t)):
    zshift[i] = t[i] * 100.0

    solution = init * np.exp(-1.0 * kz**2 * t[i] + 1j * kz * zshift[i])
    analytic = np.fft.fft(solution)[::upfactor]
    maxes[i] = np.max(analytic)

maxes[90] = np.nan
ax1.plot((zshift + 0.1) % 1.0 - 0.1, maxes, ls=":", c="k")

ax1.set_ylabel("$T$")
ax1.set_xlabel("$z$")

d = np.genfromtxt("cn_converge.dat", names=True)

mask = d["nz"] == 128

ax2.scatter(d["dt"][mask],d["l1_end"][mask], label="RK-CN2")
ax2.plot(d["dt"][mask],(d["dt"][mask]/d["dt"][mask][0])**2*d["l1_end"][mask][0], ls="--", c="r", label="$\propto \Delta t^{2}$")

d = np.genfromtxt("abbdf_converge.dat", names=True)

mask = d["nz"] == 128

ax2.scatter(d["dt"][mask],d["l1_end"][mask], label="AB-BDF3")
ax2.plot(d["dt"][mask],(d["dt"][mask]/d["dt"][mask][0])**3*d["l1_end"][mask][0], ls="--", c="k", label="$\propto \Delta t^{3}$")

ax2.legend()
ax2.set_xscale("log")
ax2.set_yscale("log")

ax2.set_xlabel("$\Delta t$")
ax2.set_ylabel("$\ell_1$")

plt.tight_layout()
plt.savefig("temporal_converge.pdf")
plt.show()
