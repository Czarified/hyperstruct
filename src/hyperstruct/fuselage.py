"""Fuselage Module.

The fuselage module operates as a stand-alone program or in conjuction with other modules.

This file contains all global variables, classes, and functions related to fuselage weight synthesis.
"""

from dataclasses import dataclass
from typing import Any
from typing import Tuple

import numpy as np

from hyperstruct import Component


@dataclass
class Cover(Component):
    """Fuselage Cover component.

    The Fuselage Cover components (aka Skins) are evalutated
    for sufficient thickness to satisfy strength, flutter,
    acoustic fatigue, and minimum thickness. Should the shell section
    be pressurized or contain fuel, they are evaluated for pressure requirements.
    The cover sizing procedure starts at minimum thicknesses and proceeds
    through a systematic check for the different criteria.
    """

    milled: bool
    """if panel is to be milled for different field vs land thicknesses."""

    L: float
    """frame spacing."""

    D: float
    """stringer spacing or panel size."""

    R: float
    """distance from the CG to the synthesis cut."""

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
    def q(self) -> float:
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
        root = float(np.sqrt(1 - self.material.nu**2))
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
            return float(7.5 * (z / 2) ** 0.113)
        elif 10 < z:
            return float(9 * (z / 10) ** 0.522)
        else:
            raise NotImplementedError("I don't know how you got here...")

    def field_thickness_block_shear(self) -> float:
        """Field thickness based on shear strength.

        Evaluate the min thickness required to satisfy block shear strength.
        """
        if self.milled:
            t_c = self.q / self.material.F_su
        else:
            t_c = self.q / (self.c_r * self.material.F_su)

        return float(t_c)

    def field_thickness_postbuckled(self, alpha: float = 45) -> float:
        """Field thickness based on critical shear flow.

        Evaluate the min thickness required to satisfy post-buckled strength.
        Postbuckled strength assumes sizing covers with diagonal tension. The
        diagonal tension angle is unknown because intermediate frame and
        stringer sizing is not known. An initial estimate of 45 degrees is used.

        Args:
            alpha: Diagonal tension angle (assumed 45 degrees)

        Returns:
            A float of field thickness.
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

            return float(t_c)

    def land_thickness_net_section(self) -> float:
        """Land thickness based on net section allowable.

        On milled panels the land thickness is checked against the net section shear allowable.
        On unmilled panels, the land thickness is simply equivalent to the field thickness.
        """
        if not self.milled:
            return float(self.t_c)
        else:
            return float(self.q / (self.c_r * self.material.F_su))

    def thickness_pressure(self, F_allow: Any = None) -> Tuple[float, float]:
        """Thicknesses based on cover pressure.

        A required thickness is evaluated to resist hoop stress,
        and then bending diaphram stress via strip theory. The
        minimum of these values is returned as the required land thickness.

        Cabin pressurization is a cyclic occurrence that subjects the cover
        to possible fatigue failure. The maximum allowable stress to prevent
        fatigue is represented as a fraction of the material ultimate
        tensile strength. The pre=programmed allowable represents a change
        from zero to peak pressure 20,000 times during the vehicle's useful
        life, with a stress concentration factor of 4.0.

        Args:
            F_allow: Allowable stress (25% yield, if not provided)

        Returns:
            A tuple of (Land thickness, Field thickness).
        """
        b = min(self.D, self.L)
        if not F_allow:
            F_allow = self.material.F_ty / 4.0

        # TODO: Lookup the vehicle-level load factors for the CG
        Nz_plus = 6
        Nz_minus = -3

        # TODO: Lookup the vehicle-level pitch accelerations
        # Assume the pitch acceleration units are radians per second per second
        Q_dot = 6.9

        # 386.0886 [in/s2] is the gravitational acceleration
        Nz_1 = Nz_plus + Q_dot * self.R * 386.0886
        Nz_2 = Nz_minus + Q_dot * self.R * 386.0886

        # TODO: Lookup the fluid density from the vehicle
        #       Should this default to air for cabin fluid?
        # 7.01 [lbs/gal] is roughly the density of JP-8
        rho = 7.01

        # TODO: Lookup the vent space pressure (atmosphere?) based on the flight condition
        # 14.7 [psi] is the atmospheric pressure at sea level.
        P_o = 14.7

        # TODO: Lookup the tank depth. Is every skin a "tank??"
        h = b / 4

        P_1 = P_o + rho * Nz_1 * h
        P_2 = P_o + rho * Nz_2 * h

        # Simple Hoop Stress
        t_1 = P_1 * self.RC / F_allow
        t_2 = P_2 * self.RC / F_allow

        # Strip Theory Edge thickness
        t_3 = (1.646 * b * P_1**0.894 * self.material.E**0.394) / F_allow**1.288
        # Strip Theory Midspan thickness
        t_4 = self.material.E**1.984 * (1.3769 * b * P_1**2.484) / F_allow**4.467

        t_min = min(t_1, t_2, t_3, t_4)

        if self.milled:
            return (float(t_3), float(t_min))
        else:
            return (float(t_min), float(t_min))

    def panel_flutter(self, mach: float, altitude: float) -> float:
        """Evaluate baseline thickness to avoid local panel flutter.

        The baseline thickness to make local panel flutter highly improbable
        is a function of dynamic pressure, Mach number, geometry, and material.
        Curve-fits for each variable are utilized, with a baseline thickness
        value returned.

        To find the final required thickness, this method must be looped
        over for all Mach and altitudes corresponding to the flight
        envelope conditions of the aircraft.

        Args:
            mach: Mach Number
            altitude: Altitude (in thousands of feet)

        Returns:
            A float of Field Thickness.
        """
        # Dynamic pressures based on standard day atmosphere.
        # Dynamic Pressure, q, in [psf]
        # Altitudes must be measured in [ft]
        if altitude <= 20:
            q = mach**2 * (1479.757 - 52.187 * altitude + 0.619868 * altitude**2)
        elif 20 < altitude <= 70:
            q = mach**2 * (
                1465.175
                - 50.76695 * altitude
                + 0.6434412 * altitude**2
                - 0.002907194 * altitude**3
            )
        elif altitude > 70:
            q = mach**2 * (199.659 / altitude) ** 4

        # Dynamic pressure is in [psf] convert to [psi]
        q = q / 144

        # Calculate the Mach Number Effect
        # The Mach Number Effect, FM, is a 3 piece function from SWEEP.
        # This function is more conservative than the AFRL baseline curve,
        # but less conservative than the NACA curve. ADA002867, pg. 100.
        if 1.0 <= mach < 1.4:
            FM = 0.4851674 + 1.66456 * (mach - 1) ** 3
        elif 1.4 <= mach <= 2.0:
            FM = (
                0.488412
                - 0.4037203 * np.cos(np.pi * (mach - 1.4) / 0.6)
                + 0.4849271 * np.sqrt(mach**2 - 1)
            )
        elif mach > 2.0:
            FM = np.sqrt(mach**2 - 1)

        # Length/Width parameter
        LW = self.L / self.D
        if LW > 10:
            # These curves are not valid for ratios greater than 10.
            # Panels with aspect ratios greater than 10 rarely exist.
            LW = 10

        # Panel Flutter Parameter, phi_B
        phi_b = 0.5551841 - 0.1686944 * LW + 0.02169992 * LW**2 - 0.000963694 * LW**3

        t_b = (phi_b * self.L) / (FM * self.material.E / q) ** (1 / 3)

        return float(t_b)

    def acoustic_fatigue(self) -> Tuple[float, float]:
        """Thickness requirements based on acoustic fatigue.

        Assumptions are:
            1.) The overall pressure level from jet engine noise is approximately
                30db higher than the random spectrum pressure level.
            2.) The accuracy of noise intensity prediction is on the order of +/-3db.
            3.) Beam theory is sufficient for modelign the shell elements.
            4.) The design to repetitive pressure intensity from jet engine noise
                can be correlated tot he material endurance limit. Polished specimen
                s-n data, K_t=1, is a sufficiently representative strength index.
            5.) The method is based on a limited amount of testing.
            6.) The method does not consider panel aspect ratio, and is, therefore,
                conservative when aspect ratios approach 1.

        The sound pressure level used is based on a material property. This value
        provides the fatigue life of 10^9 cycles for the material. The overall
        decibel level is then increased by 30, which represents jet noise instead
        of a purely random spectrum.

        Returns:
            A tuple of Land thickness, Field thickness (t_l, t_c).
        """
        # Random distribution acoustic level, that provides a material fatigue
        # life of 10^9 cycles.
        db_r = self.material.db_r
        db_oa = db_r + 30
        P = 2.9e-9 * 10 ** (db_oa / 20)
        # Note: K_c is directly hardcoded to 2.62, per Fig. 20 of ADA002867
        if self.milled:
            t_l = 2.62 * self.D * np.sqrt(P / self.material.F_en)
            t_c = 0.6 * t_l
        else:
            t_l = 2.62 * self.D * np.sqrt(P / self.material.F_en)
            t_c = t_l

        # Curvature correction
        x_l = self.D**2 / (t_l * self.RC)
        x_c = self.D**2 / (t_c * self.RC)
        # Ratio of curved panel thickness over flat panel
        f_l = 1.0794 + 0.000143 * x_l - 0.076475 * (1 / x_l) - 0.29969 * np.log(x_l)
        f_c = 1.0794 + 0.000143 * x_c - 0.076475 * (1 / x_c) - 0.29969 * np.log(x_c)

        return (float(f_l * t_l), float(f_c * t_c))


@dataclass
class MinorFrame(Component):
    """Fuselage Minor Frame Component.

    Minor Frames form part of the basic structural grid work that resists
    vehicle shear and bending loads. As opposed to Major Frames, Minor
    Frames do not redistribute concentrated loads. Minor Frames are sized
    for general shell stability, acoustic fatigue, and forced crippling
    due to postbuckled cover (skin) design.

    Minor Frames have an assumed construction, shown here? The geometry
    variables, frame depth (c), and cap flange width (b), are user input
    parameters. Thickness is the only parameter determined by sizing.

    For simplicity, the flanges are assumed to have a thickness twice
    that of the web.
    """

    c: float
    """frame depth."""

    b: float
    """cap flange width."""

    t_r: float = 0.0
    """cap flange thickness."""

    @property
    def t_w(self) -> float:
        """Web thickness."""
        return self.t_r / 2

    @property
    def area(self) -> float:
        """Cross-sectional area."""
        return 4 * self.b * self.t_r + self.c * self.t_w

    @property
    def i_xx(self) -> float:
        """Second moment of area (bending moment of inertia).

        Assumes the simplified case that t_r = 2*t_w.
        """
        return self.t_r * (
            self.b * self.c**2 + 2 * self.b**3 / 3 - self.b**2 * self.c + self.c**3 / 24
        )

    @property
    def rho(self) -> float:
        """Radius of gyration."""
        return float(
            (
                (
                    self.b * self.c**2
                    + 2 * self.b**3 / 3
                    - self.b**2 * self.c
                    + self.c**3 / 24
                )
                / (4 * self.b + self.c / 2)
            )
            ** 0.5
        )

    def general_stability(self, L: float, D: float, M: float) -> float:
        """Thickness to avoid general instability.

        The thickness that provides frame stiffness sufficient to prevent
        general instability failure is solved via the Shanley equation.

        Args:
            L: Frame Spacing
            D: Fuselage Diameter
            M: Bending moment at the cut

        Returns:
            A float of Flange thickness.
        """
        c_f = 1 / 16000
        numerator = c_f * M * D**2
        denominator = (
            self.material.E_c
            * L
            * (self.b * self.c**2 + 2 * self.b**3 * self.c / 3 + self.c**3 / 24)
        )

        return float(numerator / denominator)

    def acoustic_fatigue(self, b: float) -> float:
        """Thickness requirements based on acoustic fatigue.

        Assumptions are:
            1.) The overall pressure level from jet engine noise is approximately
                30db higher than the random spectrum pressure level.
            2.) The accuracy of noise intensity prediction is on the order of +/-3db.
            3.) Beam theory is sufficient for modelign the shell elements.
            4.) The design to repetitive pressure intensity from jet engine noise
                can be correlated tot he material endurance limit. Polished specimen
                s-n data, K_t=1, is a sufficiently representative strength index.
            5.) The method is based on a limited amount of testing.
            6.) The method does not consider panel aspect ratio, and is, therefore,
                conservative when aspect ratios approach 1.

        The sound pressure level used is based on a material property. This value
        provides the fatigue life of 10^9 cycles for the material. The overall
        decibel level is then increased by 30, which represents jet noise instead
        of a purely random spectrum.

        Args:
            b: Support spacing (frame spacing)

        Returns:
            A float of Flange thickness.
        """
        # Random distribution acoustic level, that provides a material fatigue
        # life of 10^9 cycles.
        db_r = self.material.db_r
        db_oa = db_r + 30
        P = 2.9e-9 * 10 ** (db_oa / 20)
        # Note: K_c is directly hardcoded to 7.025, per Fig. 20 of ADA002867
        t_r = 7.025 * b**0.5 * P**2 / self.material.F_en

        return float(t_r)

    def forced_crippling(self) -> float:
        """Thickness from force crippling."""
        return 0.0
