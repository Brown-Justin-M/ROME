import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.style.use("paper")

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(5.5, 5.5))

d = nc.Dataset("data.nc", "r")



pc = ax1.pcolormesh(d["x"], d["z"], d["TS"][0,:,0,:])
pc.set_rasterized(True)

ax1.set_xlabel("$x$")
ax1.set_ylabel("$z$")

pc = ax2.pcolormesh(d["x"], d["z"], d["TS"][-1,:,0,:])
pc.set_rasterized(True)

ax2.set_xlabel("$x$")
ax2.set_ylabel("$z$")

d = np.genfromtxt("converge.dat",names=True)

# for mode in np.unique(d["x_modes"]):
#     mask = d["x_modes"] == mode
#     ax3.plot(d["aspect"][mask],d["l1_end"][mask], label="$N_x=%i$" % mode)

# ax3.set_xlabel("${\Gamma_z}/{\Gamma_x}$")
# ax4.set_ylabel("$\ell_1$")
# ax3.set_xlim(-1.0, 2.1)
# ax3.set_yscale("log")
# ax3.legend()
# # ax3.legend(bbox_to_anchor=(0, -0.1), loc='upper left', borderaxespad=0.)

for aspect in np.unique(d["aspect"]):
    if aspect not in [0.5, 0.6, 0.7, 1.0, 1.5, 2.0]: continue
    mask = d["aspect"] == aspect
    ax3.plot(d["x_modes"][mask],d["l1_end"][mask], label=r"${\Gamma_z}/{\Gamma_x}=%.1f$" % aspect)
    ax4.plot(d["x_modes"][mask],d["l1_end"][mask], label=r"${\Gamma_z}/{\Gamma_x}=%.1f$" % aspect)
    # if aspect == 0.5:
    #     ax4.plot(d["x_modes"][mask],d["l1_start"][mask], ls="--", c="k", label="$N_x=%i$" % mode)

ax3.set_xlabel("$N_x$")
ax3.set_ylabel("$\ell_1$")
ax3.set_xscale("log")
ax3.set_yscale("log")
# ax3.legend()

rect = patches.Rectangle((28, -28), 580, 56, color="w", zorder=4)
ax4.add_patch(rect)
ax4.legend(loc="upper left")
ax4.set_axis_off()

plt.tight_layout()
plt.savefig("remap_converge.pdf", dpi=300)
plt.show()
