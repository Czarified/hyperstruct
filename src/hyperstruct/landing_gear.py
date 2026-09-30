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
    lcid: float
    vert: float
    drag: float
    side: float


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
