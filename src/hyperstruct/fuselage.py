"""Fuselage Module.

The fuselage module operates as a stand-alone program or in conjuction with other modules.

This file contains all global variables, classes, and functions related to fuselage weight synthesis.
"""

from dataclasses import dataclass
from dataclasses import field

import numpy as np

from hyperstruct import Material


@dataclass
class Cover:
    """Fuselage Cover component.

    The Fuselage Cover components (aka Skins) are evalutated
    for sufficient thickness to satisfy strength, flutter,
    acoustic fatigue, and minium thickness. Should the shell section
    be pressurized or contain fuel, they are evaluated for pressure requirements.
    The cover sizing procedure starts at minimum thicknesses and proceeds
    through a systematic check fo the different criteria.
    """

    material: Material = field(
        metadata={"description": "material the cover is made of"}
    )
    milled: bool = field(
        default=False,
        metadata={
            "description": "if panel is to be milled for different field vs land thicknesses."
        },
    )
    L: float = field(default=0, metadata={"description": "frame spacing"})
    D: float = field(
        default=0, metadata={"description": "stringer spacing or panel size"}
    )
    t_l: float = field(
        default=0, metadata={"description": "land, or edgeband, thickness"}
    )
    t_c: float = field(
        default=0, metadata={"description": "field, or acreage, thickness"}
    )
    RC: float = field(default=0, metadata={"description": "radius of curvature"})
    V: float = field(
        default=0, metadata={"description": "total vertical shear at the cut"}
    )
    Q: float = field(
        default=0, metadata={"description": "area moment of the bending elements"}
    )
    I: float = field(
        default=0,
        metadata={"description": "area moment of inertia of the bending elements"},
    )

    @property
    def c_r(self) -> float:
        """Rivet Factor.

        Assuming a design practice of rivet spacing four diameters apart,
        degradation due to this hole spacing is represented by this rivet factor.

        For milled panels the fastener hole degradation is assumed to be accounted for,
        and this factor is set to 1.0.
        """
        return 1.0 if self.milled else 0.75

    @property
    def q(self):
        """Evaluate the cover shear flow.

        Assuming the section masses are concentrated at the longerons
        or stringers and that the areas are equal, the equation for shear
        flow may be reduced to an equation dependent only on element location.
        """
        return self.V * self.Q / self.I

    @property
    def zee(self) -> float:
        """Panel buckling geometry parameter for shear buckling coefficient."""
        b = min(self.D, self.L)
        frac = b**2 / (self.RC * self.t_c)
        # Note, ADA002867 uses lowercase greek mu, but confirms that it is Poisson's Ratio.
        root = np.sqrt(1 - self.material.nu**2)
        return frac * root

    @property
    def k_s(self) -> float:
        """Shear buckling coefficient.

        Reference curves are displayed in ADA002867, Fig. 12.
        Coefficients are parameterized here based on curvature and geometry parameter zee.
        """
        z = self.zee
        if (self.RC > 1000) or (self.RC == 0):
            # Panel has no curvature. Flat panels are assumed to be a 7.5
            return 7.5
        elif z < 2:
            return 7.5
        elif 2 < z < 10:
            return 7.5 * (z / 2) ** 0.113
        elif 10 < z:
            return 9 * (z / 10) ** 0.522

    def field_thickness_block_shear(self) -> float:
        """Field thickness based on shear flow.

        Evaluate the min thickness required to satisfy block shear strength.
        """
        if self.milled:
            t_c = self.q / self.material.F_su
        else:
            t_c = self.q / (self.C_r * self.material.F_su)

        return t_c
