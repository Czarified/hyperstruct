"""Fuselage Module.

The fuselage module operates as a stand-alone program or in conjuction with other modules.

This file contains all global variables, classes, and functions related to fuselage weight synthesis.
"""

from dataclasses import dataclass

import numpy as np

from hyperstruct import Material


@dataclass
class Cover:
    """Fuselage Cover component.

    The Fuselage Cover components (aka Skins) are evalutated
    for sufficient thickness to satisfy strength, flutter,
    acoustic fatigue, and minimum thickness. Should the shell section
    be pressurized or contain fuel, they are evaluated for pressure requirements.
    The cover sizing procedure starts at minimum thicknesses and proceeds
    through a systematic check for the different criteria.
    """

    material: Material
    """material the cover is made of."""

    milled: bool
    """if panel is to be milled for different field vs land thicknesses."""

    L: float
    """frame spacing."""

    D: float
    """stringer spacing or panel size."""

    t_l: float = 0
    """land, or edgeband, thickness."""

    t_c: float = 0
    """field, or acreage, thickness."""

    RC: float = 0
    """radius of curvature."""

    V: float = 0
    """total vertical shear at the cut."""

    Q: float = 0
    """area moment of the bending elements."""

    I: float = 0
    """area moment of inertia of the bending elements."""

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
        """Field thickness based on shear strength.

        Evaluate the min thickness required to satisfy block shear strength.
        """
        if self.milled:
            t_c = self.q / self.material.F_su
        else:
            t_c = self.q / (self.C_r * self.material.F_su)

        return t_c

    def field_thickness_postbuckled(self, alpha: float = 45) -> float:
        """Field thickness based on critical shear flow.

        Evaluate the min thickness required to satisfy post-buckled strength.
        Postbuckled strength assumes sizing covers with diagonal tension. The
        diagonal tension angle is unknown because intermediate frame and
        stringer sizing is not known. An initial estimate of 45 degrees is used.
        """
        F_scr = (
            self.k_s
            * np.pi**2
            * self.material.E
            / (12 * (1 - self.material.nu**2))
            * (self.t_c / min(self.D, self.L)) ** 2
        )
        if F_scr > self.material.F_su:
            # Skin is not critical in stability.
            return self.t_c
        else:
            # Skin must be sized for postbuckled strength.
            K_1 = (
                self.k_s
                * np.pi**2
                * self.material.E
                / (12 * (1 - self.material.nu) ** 2 * min(self.D, self.L))
            )
            K_2 = (
                self.c_r
                * self.material.F_tu
                * np.sin(np.radians(alpha))
                * np.cos(np.radians(alpha))
                / (1 - np.sin(np.radians(alpha)) * np.cos(np.radians(alpha)))
            )
            K_3 = self.q / (1 - np.sin(np.radians(alpha)) * np.cos(np.radians(alpha)))

            t_c = (
                K_3 / (2 * K_1)
                + np.sqrt((K_3 / (2 * K_1)) ** 2 + (K_2 / (3 * K_1)) ** 3)
            ) ** (1 / 3) - (
                np.sqrt((K_3 / (2 * K_1)) ** 2 + (K_2 / (3 * K_1)) ** 3)
                - (K_3 / (2 * K_1))
            ) ** (
                1 / 3
            )

            return t_c

    def land_thickness_net_section(self) -> float:
        """Land thickness based on net section allowable.

        On milled panels the land thickness is checked against the net section shear allowable.
        On unmilled panels, the land thickness is simply equivalent to the field thickness.
        """
        if not self.milled:
            return self.t_c
        else:
            return self.q / (self.c_r * self.material.F_su)
