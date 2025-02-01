"""This is a testing sandbox script for Diagonal Tension.

This script builds a couple examples from Bruhn, and executes the
hyperstruct diagonal tension class methods to compare results.

Reference Bruhn "Analysis and Design of Flight Vehicle Structures", 1973
"""

from hyperstruct import Material
from hyperstruct.fuselage import ForcedCrippling


# Example C11.34
al_24st = Material(
    rho=0.1,
    E=10.6e6,
    E_c=10.6e6,
    nu=0.3,
    F_tu=64e3,
    F_ty=42.1e3,
    F_cy=48.3e3,
    F_su=41.0e3,
    F_bru=10.04e3,
    F_bry=89.0e3,
    F_en=20.0e3,
    db_r=116,
)

example_34 = ForcedCrippling(
    material=al_24st, c=3.0, b=1.0, construction="stringer", t_r=0.040
)

stuff = example_34.forced_crippling(
    L=15.0,
    D=65.8,
    M=1.475e6,
    Z=32.9,
    sum_z_sq=120,
    t_c=0.032,
    RC=32.9,
    f_s=6000,
    f_scr=815,
    cover_e=10.6e6,
    long_e=10.6e6,
    frame_e=10.6e6,
)

print(stuff)
