"""This script builds components for an example Medium Transport aircraft.

All the values in this file are taken from publicly available data
on the Lockheed C-130J-30. However, instead of matching the fuselage shape
exactly, a rounded rectangle will be used.

Assumptions are made for material properties, and material data taken
from either SWEEP documentation or Matweb.

Unit system is lbf, in, s

References:
https://en.wikipedia.org/wiki/Lockheed_C-130_Hercules
https://man.fas.org/dod-101/sys/ac/c-130.htm
https://www.lockheedmartin.com/content/dam/lockheed-martin/aero/documents/sustainment/csc/service-news/sn-mag-v1-v10/V2N1.pdf
"""

import matplotlib.pyplot as plt

from hyperstruct import Material
from hyperstruct import Station


# from hyperstruct.fuselage import Cover
# from hyperstruct.fuselage import ForcedCrippling
# from hyperstruct.fuselage import Fuselage
# from hyperstruct.fuselage import MajorFrame


# Some global variables for reference
FUSELAGE_LENGTH = 12 * (97 + 15) + 9
# Ignore the floor break transitions, since we're not matching the exact shape
# we'll turn this into a rounded rectangle of equivalent length and width.
FUSELAGE_DIAMETER = 12 * 14.17
CENTER_FUSE_Z_REF = 12 * 7
# A dict of FS locations where Major Frames exist and shapes change
_w, _d, _z, _r = [
    FUSELAGE_DIAMETER,
    FUSELAGE_DIAMETER,
    CENTER_FUSE_Z_REF,
    0.45 * FUSELAGE_DIAMETER,
]
FS_DICT = {
    93: [
        "Nose Station, and NLG Bay Start",
        12 * 2 + 10,
        12 * 3 + 2,
        0.5 * CENTER_FUSE_Z_REF,
        13,
    ],
    165: ["FWD Fuse, and NLG Bay End", _w, _d, _z, _r],
    245: ["FWD Fuse Mate", _w, _d, _z, _r],
    497: ["MLG Bay Start, Wing Fwd Spar", _w, _d, _z, _r],
    620: ["MLG Bay End, Wing Aft Spar", _w, _d, _z, _r],
    737: ["Empennage Mate", _w, 0.9 * _d, _z + 0.1 * _d, _r],
    941: ["Ramp Cutout End", _w, 0.5 * _d, _z + 0.5 * _d, 0.5 * _r],
    1071: [
        "Vert and Horz Stabilizer Fwd Spar",
        0.7 * _w,
        0.3 * _d,
        _z + 0.7 * _d,
        0.3 * _r,
    ],
    1200: ["Tail Station", 0.5 * _w, 0.2 * _d, _z + 0.7 * _d, 0.2 * _r],
}
#

# Build the stations geometry
stations = []
for fs, values in FS_DICT.items():
    name, width, depth, z_ref, radius = values
    stations.append(
        Station(
            orientation="FS",
            name=name,
            number=fs,
            width=width,
            depth=depth,
            vertical_centroid=z_ref,
            radius=radius,
        )
    )

# Plot the geometry
fig, axs = plt.subplots(
    nrows=3,
    ncols=3,
    figsize=(9, 9),
    layout="constrained",
    gridspec_kw=dict(hspace=0, wspace=0),
    sharex="col",
    sharey="row",
    subplot_kw=dict(box_aspect=1),
)
axs = axs.flatten()
square_side = 1.1 * _d + _z
xlim = (-square_side / 2, square_side / 2)
ylim = (0, square_side)
# print(f"x-axis width = {xlim[1]-xlim[0]:.2f}")
# print(f"y-axis depth = {ylim[1]-ylim[0]:.2f}")

for i, station in enumerate(stations):
    station.show(display=False, ax=axs[i], xlim=xlim, ylim=ylim)

for ax in fig.get_axes():
    ax.label_outer()

fig.suptitle(
    "Station Diagrams for Generic Medium Transport", fontfamily="serif", fontsize=16
)
plt.show()


# A basic aluminum material
# 2024-T3, Sheet, A-basis
al2024 = Material(
    rho=0.1,
    E=10.5e6,
    E_c=10.6e6,
    nu=0.33,
    F_tu=64e3,
    F_ty=42.1e3,
    F_cy=48.3e3,
    F_su=41.0e3,
    F_bru=10.04e3,
    F_bry=89.0e3,
    F_en=20.0e3,
    db_r=116,
)
