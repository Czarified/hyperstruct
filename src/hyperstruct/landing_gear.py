"""Landing Gear module. [ADA002858_Vol5.2].

This module analyzes and sizes landing gear for weight predictions.


Working Notes as I'm reading the reference material:
    - SWEEP has 1 main and 5 subroutines for this. Broken down into:
        -- Read/process input data, main module [LANDGR]
        -- Drag/Side/Vert loads on wheels [LGEAR]
        -- axial/normal loads on strut [LOADS]
        -- calculate weight based on loads [LGWT]
        -- bending/torsion modulus of rupture [BMOR]
        -- 3-pt interpolation for final results [LG3P]
    - Landing Gear analysis needs design parameter inputs:
         1. Takeoff/Landing weights
         2. Takeoff/Landing CGs
         3. Height of CG to ground plane
         4. Wing area
         5. Takeoff/Landing speeds
         6. Takeoff/Landing sink rates
         7. Load factor (input or calculated)
         8. Coefficient of lift at takeoff/landing
         9. Materials
        10. Gear locations in FS/BL/WL (see 3)
        11. Length/Stroke of struts
        12. Strut piston diameter (input or calcaulated)
        13. Wheel eccentricity
        14. Number of wheels
        15. Strut angles fore/aft of main/nose
        16. Tire dimensions
    - Loads are based on MIL-A-008862, and use 8 load conditions
      (2-pt, spinup, springback, braked roll, drift, unsymmetric braking,
      towing, and turning).
    - Components for sizing are broken down into:
        1. Outer Cylinder
        2. Inner Cylinder (piston)
        3. Axle
        4. Bogie
        5. Drag and Side Struts
        6. Oil
        7. Tires, Tubes, and Wheels
        8. Brakes
    - Method descriptions are documented in the order the program uses them!
    -
"""

# import logging
# from collections import namedtuple
# from copy import copy
from dataclasses import dataclass

# from typing import Dict
# from typing import Any
# from typing import NamedTuple
from typing import Tuple

# import matplotlib.pyplot as plt
import numpy as np

# from hyperstruct import Component

# from dataclasses import field


# import pandas as pd
# from matplotlib.axes import Axes
# from matplotlib.figure import Figure
# from matplotlib.lines import Line2D
# from matplotlib.patches import FancyArrowPatch
# from numpy.typing import ArrayLike
# from rich.logging import RichHandler
# from scipy.optimize import minimize_scalar


# from hyperstruct import LoadCase
# from hyperstruct import Material


#
# Classes
#


@dataclass
class GroundLoads:
    """Ground Loads base class."""

    name: str
    lcid: int
    main: dict
    nose: dict


#
# Module Functions
#


def landing_speed(
    grwt_to: float,
    grwt_l: float,
    dwt: float,
    s_w: float,
    clift_to: float,
    clift_l: float,
) -> Tuple[float, float]:
    """Calculate landing speeds at takeoff (to), and landing (l).

    Landing speeds are a function of the aircraft weight and lift. Speeds are
    calculated at 2 different weight conditions representing an aborted takeoff
    and a regular landing.

    Args:
        grwt_to (float): Gross weight at takeoff (ft/s)
        grwt_l (float): Gross weight at landing (ft/s)
        dwt (float): Aborted takeoff delta weight (lbs)
        s_w (float): Wing area (ft2)
        clift_to (float): coefficient of lift at takeoff weight
        clift_l (float): coeficcient of lift at landing weight

    Returns:
        Tuple[float, float]: landing speed at takeoff and landing weights (respectively).
    """
    vl_to = 34.776 * ((grwt_to - dwt) / (s_w * clift_to)) ** 0.5
    vl_l = 34.776 * (grwt_l / (s_w * clift_l)) ** 0.5

    return (vl_to, vl_l)


def load_factors(
    fea: float,
    ss_to: float,
    ss_l: float,
    clift_w: float,
    stroke_to: float,
    stroke_l: float,
    od_m: float,
    g: float = 32.172,
) -> Tuple[float, float]:
    """Calculate the Load Factors.

    Load Factors are calculated from the strokes, sink speeds, wing lift
    coeffciient, and the tire diameter.

    Args:
        fea (float): fraction of energy absorbed by strut.
        ss_to (float): sink speed at takeoff weight.
        ss_l (float): sink speed at landing weight.
        clift_w (float): wing lift coefficient.
        stroke_to (float): effective stroke of MLG at takeoff weight.
        stroke_l (float): effective stroke of MLG at landing weight.
        od_m (float): outer diameter of MLG tires.
        g (float, optional): gravitational constant. Defaults to 32.172.

    Returns:
        Tuple[float, float]: Load Factors at takeoff weight and landing weight.
    """
    ng_to = (
        (1 - fea)
        * (ss_to**2 / (2 * g) + (1 - clift_w) * (0.98 * stroke_to + 0.08 * od_m / 12))
        / (0.8 * stroke_to)
    ) + clift_w

    ng_l = (
        (1 - fea)
        * (ss_l**2 / (2 * g) + (1 - clift_w) * (0.98 * stroke_l + 0.08 * od_m / 12))
        / (0.8 * stroke_l)
    ) + clift_w

    return (ng_to, ng_l)


def piston_diameters(
    grwt_to: float, cg_to: float, fs_n: float, fs_m: float, strut_m: int
) -> Tuple[float, float]:
    """MLG and NLG Piston Diameters.

    The Main Gear piston diameter is a function of static load.
    The Nose Gear piston diameter is simply a ratio of the main.
    A different function is used for static loads over 77,295 lbs.

    Args:
        grwt_to (float): gross weight at takeoff.
        cg_to (float): center of gravity at takeoff.
        fs_n (float): fuselage station of nose gear.
        fs_m (float): fuselage station of main gear.
        strut_m (int): number of main gear struts.

    Returns:
        Tuple[float, float]: Main and Nose gear piston diameters.
    """
    sw = grwt_to * np.abs((cg_to - fs_n) / (fs_m - fs_n)) / strut_m

    if sw > 77295:
        dp_m = ((4 * sw) / (15000 * np.pi)) ** 0.5
    else:
        if sw < 5542:
            acm = 187.5
            bcm = 380.0
        elif sw < 33819:
            acm = 126.7
            bcm = 545.0
        elif sw <= 77295:
            acm = 95.6
            bcm = 720.0

        aa = -0.333 * (bcm / acm) ** 2
        bb = 2 / 27 * (bcm / acm) ** 3 - (4 * sw) / (np.pi * acm)
        radpd = (bb**2 / 4 + aa**3 / 27) ** 0.5

        dp_m = (-bb / 2 + radpd) ** 0.333 + (-bb / 2 - radpd) ** 0.333 - bcm / (3 * acm)

    dp_n = 0.6 * dp_m

    return (dp_m, dp_n)


def wtt_weights(
    od_m: float, od_n: float, w_m: float, w_n: float, ws_m: float, ws_n: float
) -> Tuple[float, float]:
    """Wheel, Tire, and Tube weights.

    The wheel, tire, and tube weights are calculated from the
    width and diameter of the wheels. These are statistical
    methods. The original data source is unknown.

    45% of the wheel, tire, and tube weight is in the wheels;
    therefore the total wheel, tire, and tube weights can be computed.

    Args:
        od_m (float): outer diameter of the main tires
        od_n (float): outer diameter of the nose tire(s)
        w_m (float): width of the main tires
        w_n (float): width of the nose tire(s)
        ws_m (float): wheels per strut on main gear
        ws_n (float): wheels per strut on nose gear

    Returns:
        Tuple[Tuple[float, float], Tuple[float, float]]: Main and Nose gear tuples
        containing the wheel and tube/tire weights, respectively.
        `((MLG_wheel, MLG_tubetire), (NLG_wheel, NLG_tubetire))`
    """
    # Weight per wheel of main gear (wheel, tire, and tube)
    wtt_m = 0.425 * od_m * w_m + 0.00023 * ((od_m * w_m) / 100) ** 7
    # Weight per wheel of nose gear (wheel, tire, and tube)
    wtt_n = 0.4 * od_n * w_n + 0.00024 * ((od_n * w_n) / 100) ** 8

    # Assume 45% is the wheel
    # The weight *per aircraft* of main and nose wheels
    wheel_m = 0.45 * ws_m * wtt_m**2
    wheel_n = 0.45 * ws_n * wtt_n

    # The weight *per aircraft* of main and nose tube+tires
    tt_m = 1.222 * wheel_m
    tt_n = 1.222 * wheel_n

    return ((wheel_m, tt_m), (wheel_n, tt_n))


def brake_weight(grwt_to: float, vl_to: float) -> float:
    """Brake Weights.

    The weight of brakes per aircraft is calculated via statistical
    correlation to the takeoff weight, and landing speed.

    Args:
        grwt_to (float): gross weight at takeoff
        vl_to (float): landing speed at takeoff weight

    Returns:
        float: total weight of brakes on the aircraft
    """
    brakes = 0.010783 * grwt_to * vl_to**2 * 0.00000408
    return brakes


def rotating_inertia(
    od_m: float, tt_m: float, w_m: float, brakes: float, wheel_m: float, strut_m: float
) -> float:
    """Polar mass moment of inertia for main gear.

    The inertia for the main gear wheels, tires, tubes, and brakes is calculated
    from the wheel, tire, tube, and brake weights and the tire dimensions.

    Args:
        od_m (float): outer diameter of the main gear tires
        tt_m (float): weight per aircraft of main gear tube/tires
        w_m (float): width of main gear tires
        brakes (float): total weight of the brakes on the aircraft
        wheel_m (float): weight per aircraft of main gear wheels
        strut_m (float): number of main gear struts

    Returns:
        float: inertia per strut of main gear (slug-ft2)
    """
    g = 32.172
    iw_m = (
        (od_m / (12 * 2.52)) ** 2 * tt_m
        + ((od_m - 1.818 * w_m) / (12 * 2.5)) ** 2 * (0.65 * brakes + wheel_m)
    ) / (strut_m * g)

    return iw_m


def strut_loads(
    vf: float,
    df: float,
    sf: float,
    theta_1: float,
    theta_2: float,
    is_main: bool = True,
) -> Tuple[float, float, float]:
    """Axial and Normal strut loads.

    This function calculates the axial and normal strut loads
    based on the ground reactions at the wheels, and strut
    angles. The normal load is the resultant shear load on the strut.

    Args:
        vf (float): vertical force from the ground reaction
        df (float): drag force from the ground reaction
        sf (float): side force from the ground reaction
        theta_1 (float): fore-aft angle of strut, radians
        theta_2 (float): lateral angle of strut, radians
        is_main (bool): if the gear in question is main or nose. Defaults to True (main gear).

    Returns:
        Tuple[float, float, float]: resultant load, axial load, and normal load
    """
    rload = np.sqrt(vf**2 + df**2 + sf**2)

    # Direction cosines of the resultant load
    crv = vf / rload
    crfa = df / rload
    crl = sf / rload

    if is_main:
        # Direction cosines of the main gear struts.
        # Cosine of angle between strut and vertical
        csv = np.cos(
            np.arctan(np.cos(theta_1) ** (-2) + np.cos(theta_2) ** (-2) - 2) ** 0.5
        )
        # Cosine of angle between strut and fore-aft
        csfa = np.cos(
            np.arctan(np.sin(theta_1) ** (-2) + np.cos(theta_2) ** (-2) - 2) ** 0.5
        )
        # Cosine of angle between strut and lateral
        # Is this one actually the same as csv? The manual repeats the same formula.
        csl = np.cos(
            np.arctan(np.cos(theta_1) ** (-2) + np.cos(theta_2) ** (-2) - 2) ** 0.5
        )
    else:
        # Direction cosines of the nose gear struts.
        csv = np.cos(theta_1)
        csfa = np.sin(theta_2)
        csl = 0

    # The combined angle between the resultant load and the strut.
    theta = np.arccos(csv * crv + csfa * crfa + csl * crl)

    aload = rload * np.cos(theta)
    pload = rload * np.sin(theta)

    return (rload, aload, pload)


# Landing and Ground Loads Methods
# ---------------------------------
# This section has all the unique ground loads (2-pt, spinup, springback,
# braked roll, drift, unsymmetric braking, towing, and turning).

# The ground reactions on the wheels (VF, DF, and SF) for each load condition are
# determined in accordance with the procedure outlined in MIL-A-008862. After
# the loads have been determined, the program then - except for the spring-back
# condition - uses the method described in the `strut_loads` function to find
# the axial and normal components.


def two_point_landing(
    ng_to: float,
    ng_l: float,
    cl_w: float,
    grwt_to: float,
    grwt_l: float,
    a_to: float,
    a_l: float,
    dist: float,
    dwt: float = 0.0,
) -> Tuple[GroundLoads, GroundLoads]:
    """2-PT Vert Landing.

    The vertical load on the wheels at the 2-pt landing condition is the
    maximum vertical load. The nose gear load is determined as a ratio of
    the main gear load. The landing loads are determined at both takeoff
    and landing vehicle weights. The drag load is set to one quarter of
    the vertical load. The side load is assumed to be zero.

    Args:
        ng_to (float): load factor at takeoff
        ng_l (float): load factor at landing
        cl_w (float): wing lift coefficient
        grwt_to (float): gross weight at takeoff
        grwt_l (float): gross weight at landing
        a_to (float): distance from CG to main gear, at takeoff
        a_l (float): distance from CG to main gear, at landing
        dist (float): distance from main to nose
        dwt (float): aborted takeoff delta weight. Defaults to 0.

    Returns:
        Tuple: Load collectors for takeoff and landing weights
    """
    # Maximum vertical load on main gear for takeoff and landing
    vmxmg_to = (1.5 * (ng_to - cl_w) * (grwt_to - dwt)) / 2
    vmxmg_l = (1.5 * (ng_l - cl_w) * (grwt_l)) / 2
    # Drag force
    dmxmg_to = 0.25 * vmxmg_to
    dmxmg_l = 0.25 * vmxmg_l
    # Side force
    smxmg_to = 0
    smxmg_l = 0

    # Maximum vertical loads on the nose gear
    vmxng_to = 2 * vmxmg_to * (a_to / dist)
    vmxng_l = 2 * vmxmg_l * (a_l / dist)
    # Drag force
    dmxng_to = 0.25 * vmxng_to
    dmxng_l = 0.25 * vmxng_l
    # Side force
    smxng_to = 0
    smxng_l = 0

    takeoff = GroundLoads(
        name=f"2-PT Landing, Aborted TO, W={grwt_to:.0f}lbs",
        lcid=101,
        main={"vf": vmxmg_to, "df": dmxmg_to, "sf": smxmg_to},
        nose={"vf": vmxng_to, "df": dmxng_to, "sf": smxng_to},
    )
    landing = GroundLoads(
        name=f"2-PT Landing, W={grwt_l:.0f}lbs",
        lcid=102,
        main={"vf": vmxmg_l, "df": dmxmg_l, "sf": smxmg_l},
        nose={"vf": vmxng_l, "df": dmxng_l, "sf": smxng_l},
    )

    return (takeoff, landing)
