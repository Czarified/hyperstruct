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
# from dataclasses import dataclass
# from dataclasses import field

# from typing import Dict
# from typing import Any
# from typing import NamedTuple
from typing import Tuple

# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# from matplotlib.axes import Axes
# from matplotlib.figure import Figure
# from matplotlib.lines import Line2D
# from matplotlib.patches import FancyArrowPatch
# from numpy.typing import ArrayLike
# from rich.logging import RichHandler
# from scipy.optimize import minimize_scalar

# from hyperstruct import Component
# from hyperstruct import LoadCase
# from hyperstruct import Material

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
