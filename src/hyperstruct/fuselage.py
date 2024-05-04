"""Fuselage Module.

The fuselage module operates as a stand-alone program or in conjuction with other modules.

This file contains all global variables, classes, and functions related to fuselage weight synthesis.
"""

from dataclasses import dataclass
from dataclasses import field

import numpy as np


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
    C_r: float = field(
        default=0.75,
        metadata={
            "description": "Assuming a design practice of rivet spacing four diameters apart, \
            degradation due to this hole spacing is represented by this rivet factor."
        },
    )

    @property
    def q(self):
        """Evaluate the cover shear flow.

        Assuming the section masses are concentrated at the longerons
        or stringers and that the areas are equal, the equation for shear
        flow may be reduced to an equation dependent only on element location.
        """
        return self.V * self.Q / self.I

    @property
    def zee(self):
        """Panel buckling parameter for shear buckling coefficient."""
        b = min(self.D, self.L)
        frac = b**2 / (self.RC * self.t_c)
        # Note, ADA002867 uses lowercase greek mu, but confirms that it is Poisson's Ratio.
        root = np.sqrt(1 - self.material.nu**2)
        return frac * root

    def min_thickness_block_shear(self) -> None:
        """Evaluate the min thickness requires tos atisfy block shear strength."""
        if self.milled:
            t_c = self.q / self.material.F_su
        else:
            t_c = self.q / (self.C_r * self.material.F_su)

        return t_c
