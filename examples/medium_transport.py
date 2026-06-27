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
import numpy as np

from hyperstruct import Material, LoadCase
from hyperstruct import Station
from hyperstruct import composite_cg

# from hyperstruct.fuselage import Cover
# from hyperstruct.fuselage import ForcedCrippling
from hyperstruct.fuselage import Fuselage
from hyperstruct.fuselage import MajorFrame, MinorFrame, Cover, Longeron


# from hyperstruct.fuselage import MinorFrame


# Some global variables for reference
FUSELAGE_LENGTH = 12 * (97 + 15) + 9
# Ignore the floor break transitions, since we're not matching the exact shape
# we'll turn this into a rounded rectangle of equivalent length and width.
FUSELAGE_DIAMETER = 12 * 14.17
CENTER_FUSE_Z_REF = 12 * 7 + 20
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
    station.show(display=False, axes=axs[i], xlim=xlim, ylim=ylim)

for ax in fig.get_axes():
    ax.label_outer()

fig.suptitle(
    "Station Diagrams for Generic Medium Transport", fontfamily="serif", fontsize=16
)
# plt.show()


# An initial estimate at weight distribution based on CG target
target_weight = 73000
point_weights = np.array(
    [
        # Weight , FS
        [-1200, 93],
        [-6500, 165],
        [-9000, 245],
        [-15000, 497],
        [-16000, 620],
        [-11000, 737],
        [-7600, 941],
        [-4200, 1071],
        [-2500, 1200],
    ]
)
_w, _cg = composite_cg(point_weights)
print(f"Target Weight = {target_weight:d} [lbs]")
print(f" Total Weight = {_w:d} [lbs]")
print(f"           CG = {_cg:.2f} [in]")


#
# Landing Gear Loads Calculation
#       Taxi, WC=73k, 2.0g
#
FNZ0 = 2.0
XNGG = 165
XMGG = 620
XCG = _cg
DGW = target_weight
Rmg = FNZ0 * DGW * (XCG - XNGG) / (XMGG - XNGG)
Rng = FNZ0 * DGW - Rmg

# Sink Speed is 10 ft/s, convert to in/s
# SSPD = 12 * 10.0
# Shock Strut Stroke
# STKE = 26
# g = 386     # in/s2

# Assume no Airloads for now
PZN = 0  # Forebody lift
XCPN = 170  # Forebody Center of pressure
PZWB = 0  # Wing outer panel lift
XCPW = 620  # Coord of outer wing panel
PZBW = 0  # Body lift in presence of wing
XCPB = 600  # Coord of of body lift
PZH = 0  # Horz Tail Lift
XCPH = 1200  # Coord of hTail lift

# Vehicle Pitch Inertia, simplified. See pg 56 fo ADA002867
TIYY = np.sum(point_weights[:, 0] * (XCG - point_weights[:, 1]) ** 2)
print(f"Vehicle pitch inertia = {TIYY:.3e}")

print(f"\nLoads for {FNZ0}g Taxi:")
print(25 * "-")
print(f" Rmg = {Rmg:.3e} [lbs]")
print(f" Rng = {Rng:.3e} [lbs]")
print(f"FNZO = {FNZ0:.2f} [g]")
diff = FNZ0 * DGW - Rmg - Rng
print(
    f"Balance of Vertical Forces: {FNZ0:.1f}*{DGW:.2e} - {Rmg:.2e} - {Rng:.2e} = {diff:.1f}\n"
)


#
# Turn the Stations and Loads into Frames
#

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


inertia_loads = FNZ0 * point_weights[:, 0].flatten()

gear_loads = [
    #   y,    z,     V,   H,   M
    np.array(
        [  # FS 165
            [20.0, 17.0, Rng / 2, 0.0, 0.0],  # NLG LH Mount
            [-20.0, 17.0, Rng / 2, 0.0, 0.0],  # NLG RH Mount
        ]
    ),
    np.array(
        [  # FS 620
            [FUSELAGE_DIAMETER / 2, 47.0, Rmg / 2, 0.0, 0.0],  # MLG LH Mount
            [-FUSELAGE_DIAMETER / 2, 47.0, Rmg / 2, 0.0, 0.0],  # MLG RH Mount
        ]
    ),
]

frames = {}
for station in stations:
    if station.number == 165:
        load = gear_loads[0]
    elif station.number == 620:
        load = gear_loads[1]
    else:
        load = np.zeros((2, 5))

    frames[station.number] = MajorFrame(
        material=al2024,
        fs_loc=station.number,
        loads=load,
        geom=station,
        fd=6.0,
    )

nlg_frame = frames.pop(165)
mlg_frame = frames.pop(620)
# Don't need all the frames since only 2 of them have loads
frames = (nlg_frame, mlg_frame)
fig1, ax1 = nlg_frame.show()
fig2, ax2 = mlg_frame.show()


#
# Compile the frames and stations into a fuselage
#
w_fus = np.column_stack(
    (point_weights[:, 1], FNZ0 * point_weights[:, 0], np.zeros((len(point_weights),)))
)
w_fc = np.zeros(3)
p_air = np.zeros(3)
ext_loads = [arr[:, 2].sum() for arr in gear_loads]
p_ext = np.column_stack(
    (np.array([[165], [620]]), np.transpose(ext_loads), np.zeros((2,)))
)
print("Fuselage Frame Weights:")
print(25 * "-")
print(w_fus)
print(11 * " " + f"{np.sum(w_fus[:, 1]):.2f}")
print("\nFuselage Frame Loads:")
print(25 * "-")
print(p_ext)
print(17 * " " + f"{np.sum(p_ext[:, 1]):.2f}")

cover_model = Cover(
    material=al2024, milled=False, L=30, D=20, R=1, RC=25
)
long_model = Longeron(
    material=al2024, b=2.0, t_s=0.1, k=0.8
)
frame_model = MinorFrame(
    material=al2024, c=4.0, b=3.0, construction="longeron"
)

fuse = Fuselage(
    stations=stations, 
    major_frames=frames,
    construction="longeron",
    cover_model=cover_model,
    long_model=long_model,
    frame_model=frame_model
)
loads = fuse.net_loads(w_fus, w_fc, p_air, p_ext)
fig, (ax1, ax2) = fuse.vmt_diagram(w_fus, w_fc, p_air, p_ext)
_ = fig.suptitle(f"{FNZ0:.1f}g Taxi, WC=73kip, xCG={XCG:.0f}[in]")

with np.printoptions(precision=3):
    print("     FS     , P     ,   M_ext   ,   V     ,    M_int")
    print(loads)
    print("\n\n")

x, v, m = fuse.lookup_loads(x=400, loads=loads)
print(f"FS{x}: V={v / 1000:.1f}[kip], M={m:.2e}[in-lbs]")

_ = ax1.plot(x, v, marker="^", color="k")
_ = ax2.plot(x, m, marker="^", color="k", label="Analysis Point")
_ = ax2.legend()

plt.show()


#
#   S I Z I N G
#

lc = LoadCase(fuse_loads=loads, lcid=31, name="3g Taxi", mach=0.1, altitude=0.0)

# Major Frames
print("Major Frame Sizing:")
for frame in frames:
    frame.synthesis()
    print(f"   FS {frame.fs_loc}: {frame.weight:.1f}[lbf]")

fuse.synthesis(loadcase=lc)