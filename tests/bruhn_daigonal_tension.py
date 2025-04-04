"""This is a testing sandbox script for Diagonal Tension.

This script builds a couple examples from Bruhn, and executes the
hyperstruct diagonal tension class methods to compare results.

Reference Bruhn "Analysis and Design of Flight Vehicle Structures", 1973
"""

from hyperstruct import Material
from hyperstruct.fuselage import ForcedCrippling


# Example C11.34
# Everything here is 24st material
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

print("Bruhn Example C11.34:")
example_34 = ForcedCrippling(
    d=15,
    h=7.25,
    c=3.0,
    b=1.0,
    construction="stringer",
    frame_material=al_24st,
    cover_material=al_24st,
    long_material=al_24st,
    t_r=0.040,
)

t_r, f_st, f_rg = example_34.forced_crippling(
    D=65.8, M=1.475e6, Z=32.9, sum_z_sq=120, t_c=0.032, RC=30.0, f_s=6000, f_scr=815
)

# Bruhn values:
t_r_b = 0.040
f_st_b = 26490
f_rg_b = 7610
print(f" Minimum frame cap thickness = {t_r:.3f} vs {t_r_b:.3f} [in]")
print(f"Secondary stress in stringer = {f_st:.0f} vs {f_st_b:.0f} [psi]")
print(f"    Secondary stress in ring = {f_rg:.0f}  vs  {f_rg_b:.0f} [psi]")


# Example C11.37
print("\nBruhn Example C11.37:")
al_75st = Material(
    rho=0.1,
    E=10.6e6,
    E_c=10.5e6,
    nu=0.3,
    F_tu=64e3,
    F_ty=42.1e3,
    F_cy=70.0e3,
    F_su=41.0e3,
    F_bru=10.04e3,
    F_bry=89.0e3,
    F_en=20.0e3,
    db_r=116,
)

example_37 = ForcedCrippling(
    d=5.0,
    h=41.8,
    c=1.5,
    b=0.5,
    construction="longeron",
    frame_material=al_75st,
    cover_material=al_75st,
    long_material=al_75st,
    t_r=0.032,
)

t_r, f_st, f_rg = example_34.forced_crippling(
    D=60.0,
    M=2.135e6,
    Z=39.9,
    sum_z_sq=sum([23.1**2, 23.1**2, 23.1**2, 23.1**2]),
    t_c=0.025,
    RC=30.0,
    f_s=4240,
    f_scr=313,
)
# Bruhn values:
t_r_b = 0.032
f_st_b = 37200 + 12780
f_rg_b = 27900
print(f" Minimum frame cap thickness = {t_r:.3f} vs {t_r_b:.3f} [in]")
print(f"Secondary stress in longeron = {f_st:.0f} vs {f_st_b:.0f} [psi]")
print(f"    Secondary stress in ring = {f_rg:.0f}  vs {f_rg_b:.0f} [psi]")
