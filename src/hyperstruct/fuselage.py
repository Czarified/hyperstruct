"""Fuselage Module.

The fuselage module operates as a stand-alone program or in conjuction with other modules.

This file contains all global variables, classes, and functions related to fuselage weight synthesis.
"""

from dataclasses import dataclass
from typing import Any
from typing import Tuple

import numpy as np
from scipy.optimize import minimize_scalar

from hyperstruct import Component
from hyperstruct import Station


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

    construction: str
    """Construction method. ('stringer' or 'longeron')"""

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

    def ring_allowable_stress(
        self, RC: float, K: float, t_c: float, frame_e: float, cover_e: float
    ) -> Tuple[float, float]:
        """Allowable ring stress for the Forced Crippling method.

        Args:
            RC: Radius of curvature
            K: Diagonal tension factor
            t_c: Cover thickness
            frame_e: Frame Young's Modulus
            cover_e: Cover Young's Modulus

        Returns:
            A tuple of (F_RG, G), Frame allowable stress,
            and emiprical factor from Bruhn.
        """
        # Allowable ring frame stress
        if RC <= 151.586:
            N = (
                (18695 + 75.238 * RC)
                * K ** (2 / 3)
                * (self.t_r / t_c) ** (1 / 3)
                * (frame_e / cover_e) ** (1 / 9)
            )
        elif RC > 151.586:
            N = (
                30100
                * K ** (2 / 3)
                * (self.t_r / t_c) ** (1 / 3)
                * (frame_e / cover_e) ** (1 / 9)
            )

        G = self.material.F_cy * 1088 - 5 / (
            5.88 * (self.material.F_cy / frame_e + 0.002) ** 0.5
        )

        F_RG = N * G

        return (F_RG, G)

    def diagonal_factor(
        self, D: float, L: float, t_c: float, RC: float, f_s: float, f_scr: float
    ) -> float:
        """Calculate the diagonal tension factor.

        There's two rules here. The dimension ratios are capped at 2, and
        we need to use the larger ratio depending on the convention of construction.
        Typically, for Stringer systems D > L, but for Longerons L > D.

        Args:
            D: Fuselage Diameter
            L: Frame Spacing
            t_c: Cover thickness
            RC: Side panel radius of curvature
            f_s: Shear stress in Cover at cut
            f_scr: Critical shear buckling strength of Cover

        Returns:
            Diagonal Tension factor, K, from empirical relations.

        Raises:
            ValueError: An error occurred comparing values of D and L.
        """
        if D >= L:
            if D / L > 2:
                k = np.tanh((0.5 + 300 * (t_c / RC * 2)) * np.log10(f_s / f_scr))
            else:
                k = np.tanh((0.5 + 300 * (t_c / RC * D / L)) * np.log10(f_s / f_scr))
        elif L > D:
            if L / D > 2:
                k = np.tanh((0.5 + 300 * (t_c / RC * 2)) * np.log10(f_s / f_scr))
            else:
                k = np.tanh((0.5 + 300 * (t_c / RC * L / D)) * np.log10(f_s / f_scr))
        else:
            raise ValueError(
                "Frame Spacing (L) and Fuselage Diameter (D) cannot be properly compared."
            )

        return float(k)

    def forced_crippling(
        self,
        L: float,
        D: float,
        M: float,
        Z: float,
        sum_z_sq: float,
        t_c: float,
        RC: float,
        f_s: float,
        f_scr: float,
        cover_e: float | None = None,
        long_e: float | None = None,
    ) -> float:
        """Thickness from forced crippling.

        If the covers on the fuselage structure are allowed to buckle under shear loads,
        in the buckled state these loads are supported by diagonal tension stresses.
        Axial loads produced by cover tension field are reacted by stiffening elements
        (minor frames, stringers) that bound the panel. Since shear loads are maximum at
        the midpoint of the side panel, this condition is evaluated for elements on the
        side sectors of the shell.

        Basic formulations for the prevention of forced crippling failure due to the
        postbuckled design are taken directly from Kuhn[4] and Bruhn [1]. These methods
        have been modified to account for the condition where different materials are
        selected for cover, longeron, and minor frame design.

        Cover sizing is established to satisfy strength and other criteria within the Cover
        class. Shear stress based on this thickness is compared against the critical shear stress
        to determine whether the panel is critical for postbuckled strength. At this point
        in the sizing process, longeron (stringer) area has not been established. The
        longitudinal member area is required to define panel boundary support constraints.
        A first approximation is made for longeron area base on vehicle bending loads.

        Cover sizing is established to satisfy strength and other criteria within the Cover
        class. Shear stress based on this thickness is compared against the critical shear stress
        to determine whether the panel is critical for postbuckled strength. At this point
        in the sizing process, longeron (stringer) area has not been established. The
        longitudinal member area is required to define panel boundary support constraints.
        A first approximation is made for longeron area base on vehicle bending loads.

        Args:
            L: Frame Spacing
            D: Fuselage Diameter
            M: Bending moment at the cut
            Z: coordinate of the extreme fiber (fuselage half-depth at cut)
            sum_z_sq: sum of longeron coordinates squared
            t_c: Cover thickness
            RC: Side panel radius of curvature
            f_s: Shear stress in Cover at cut
            f_scr: Critical shear buckling strength of Cover
            cover_e: Cover Young's Modulus
            long_e: Longeron Young's Modulus

        Returns:
            A bunch of floats?

        Raises:
            AttributeError: The construction attribute is not properly defined.
            ValueError: D and L cannot be successfully compared.
        """
        # Longeron/Stringer area
        A_s = M * Z / (self.material.F_cy * sum_z_sq)

        K = self.diagonal_factor(D, L, t_c, RC, f_s, f_scr)

        # Angle of the folds is given an initial approximation.
        # This is probably good enough for weights sizing, but a more exact solution involves iteration.
        # This is the radians equivalent of 45 degrees.
        alpha = np.pi / 4  # [rad]

        # The forced crippling method is highly dependant on empirical relations.
        # As such, the specific relations we use are dependant on the construction method.
        # This method must be chosen from the user, but they can use it as a conceptual
        # sizing variable permutation if they wish.
        if not cover_e:
            cover_e = self.material.E_c

        if not long_e:
            long_e = self.material.E_c

        frame_e = self.material.E_c

        # Use the same K value as before.
        alpha_ratio = K**0.25

        # A is the curve fit ordinate, Fig. 24 of ADA002867
        A = (D / RC * np.sqrt(cover_e / f_s)) / np.sqrt(
            1
            + 0.5
            * (D * t_c * cover_e / (A_s * long_e) + L * cover_e / (self.area * frame_e))
        )

        alpha_pdt = (np.pi / 4 + 0.1443 * A) / (1 + 0.175622 * A + 0.013411 * A**2)
        alpha = alpha_pdt * alpha_ratio

        # Secondary ring load: axial load in ring due to shear stress
        P_rg = f_s * K * t_c * L * np.tan(alpha)

        # Effective ring area
        A_erg = self.area / (1 + (self.c / (2 * self.rho)) ** 2)

        if self.construction == "stringer":
            ####
            # STRINGER CONSTRUCTION
            ####

            # Secondary stringer load: axial load in stringer due to shear stress
            # P_st = f_s * K * np.tan(alpha) ** (-1)

            # Effecting ring and cover area
            A_edt = A_erg + 0.5 * L * t_c * (1 - K) * (cover_e / frame_e)

            # Stress in the ring frame
            f_rg = P_rg / A_edt
            # Secondary stringer stress: axial load due to shear stress
            # f_st = P_st / A_edt

        elif self.construction == "longeron":
            ####
            # LONGERON CONSTRUCTION
            ####

            # Secondary stringer load: axial load in stringer due to shear stress
            # P_st = f_s * K * t_c * D / 2 * np.tan(alpha) ** (-1)

            A_est = A_s / (1 + (self.c / (2 * self.rho)) ** 2)

            # Effecting ring and cover area
            A_edt = A_est + 0.5 * D * t_c * (1 - K) * (cover_e / frame_e)

            # Stress in the ring frame
            f_rg = P_rg / A_edt
            # Secondary longeron stress: axial load in longeron due to shear stress
            # f_st = P_st / A_edt

        else:
            raise AttributeError(
                "The MinorFrame attribute 'construction' is not properly defined. Acceptable options are 'stringer' or 'longeron'"
            )

        # Max stress in the frame is based on an empirical relation from Bruhn.
        # This is dependent on the ratio of frame/longeron spacing.
        # This relation applies for both the applied stress, f, and the allowable stress, F.
        if L / D <= 1.2:
            frgmax_frg = 1 + 0.78 * (1 - K) - 0.65 * L / D * (1 - K)
        elif L / D > 1.2:
            frgmax_frg = 1.0
        else:
            raise ValueError(
                "Frame Spacing (L) and Fuselage Diameter (D) cannot be properly compared."
            )

        frgmax = frgmax_frg * f_rg

        F_RG, G = self.ring_allowable_stress(RC, K, t_c, frame_e, cover_e)
        H = A * G

        # Iterating for cap flange thickness, t_r
        x_c = 0.5 * (1 - K) * cover_e / self.material.E_c
        x_b = (4 * self.b + self.c / 2) / (
            (1 + (self.c / (2 * self.rho)) ** 2) * L * t_c
        )
        x_a = ((f_s * np.tan(alpha) / H) * (frgmax / F_RG)) ** 3 * t_c * K

        # An initialization for optimizing
        err = 100
        while err > 0.1:
            t_r2 = self.t_r - (
                (self.t_r * (self.t_r * x_b + x_c) ** 3 + x_a)
                / (
                    3 * x_b * self.t_r * (self.t_r * x_b + x_c) ** 2
                    + (self.t_r * x_b + x_c) ** 3
                )
            )
            err = t_r2 - self.t_r
            self.t_r = t_r2

        return float(t_r2)


@dataclass
class Longeron(Component):
    """Fuselage longitudinal member component.

    Longitudinal Members include both Stringers and Longerons. The Longitudinal Member
    sizing is determined to satisfy minimum area, forced crippling, bending strength,
    and stiffness requirements. The methods account for differences in Cover, Longeron,
    and MinorFrame materials, and the effects of cutouts.

    The flange width is set up to equal the web height, and the thickness is assumed
    to be constant across the section. This primarily simplifies
    some of the calculations, and should be sufficient for weight estimates.
    """

    b: float
    """Web Height, and Flange Width."""

    t_s: float
    """Longeron thickness."""

    k: float
    """Inner flange proportion of web height."""

    @property
    def area(self) -> float:
        """Cross-sectional area."""
        return self.t_s * self.b * (3 + self.k)

    @property
    def e(self) -> float:
        """Eccentricity."""
        return self.b * (0.5 + self.k / (3 + self.k))

    @property
    def i_xx(self) -> float:
        """Second moment of area (moment of inertia), strong-axis."""
        return (
            self.t_s
            * self.b
            * (
                2 * self.e**2
                + self.b**2 / 12
                + (0.5 * self.b - self.e) ** 2
                + self.k * (self.b - self.e) ** 2
            )
        )

    @property
    def rho(self) -> float:
        """Radius of gyration."""
        return float(np.sqrt(self.i_xx / self.area))

    @property
    def area_effective(self) -> float:
        """Effective area."""
        return self.area / (1 + (self.e / self.rho) ** 2)

    def bending_strength(
        self,
        M_ext: float,
        t: float,
        d: float,
        I_t: float,
        l_p: float,
        rtu: float,
        A_s: float,
        I_a: float,
    ) -> float:
        """Thickness required to satisfy bending strength.

        Longitudinal member sizing is dependent on the contribution of all
        copmonents that resist bending loads. This method accounts for the difference
        in behavior of cover elements under tension load versus the behavior in
        compression. Cover material, should it differ from longeron material, can
        also have different strength and elastic properties.

        The assumption that plane sections remain plane simplifies the estimating
        approach. The practice of using a Design Ultimate Factor of Safety of 1.5
        in analysis and metal deformation characteristics results in stresses at
        limit load occurring in the elastic range.

        The stress at any point on the shell is then proportional to the extreme
        fiber stress according to the relationship of vertical coordinate versus
        extreme fiber coordinate.

        Bending moment is assumed to be reactied by an internal coupled force system.
        Thus, in the case of down-bending, the uper half of the shell sustains tension
        loads, and the lower half, compression loads; half of the moment is reacted
        in each half. Covers are totally effective in the tension sector of the
        fuselage. Cutouts eliminate cover contributions; proximity to cutouts degrade
        the effectiveness of the cover. The width of cutouts at other synthesis cuts
        combined with longitudinal displacement and hsear lag slope of 2 to 1 is used
        to determine the apparent effective width.

        Args:
            M_ext: External moment at the cut
            t: cover thickness
            d: fuselage depth
            I_t: cover moment of inertia, as a function of thickness
            l_p: cover peripheral length
            rtu: apparent panel degradation due to cutout proximity
            A_s: stringer/longeron area
            I_a: side stringer moment of inertia aa a funciton of area

        Returns:
            Float of area, A_1, that satisfies the bending strength requirement
        """
        # Max allowable extreme fiber stresses
        # Longeron/Stringer
        f_max_l = 0.9 * self.material.F_cy
        # Cover
        f_max_c = 0.76 * self.material.F_tu

        # Moment reacted by the upper cover in tension
        M_c = t * (f_max_c / (d / 2)) * I_t * (l_p - rtu / l_p)

        # Moment reacted by the upper side panel stringers
        M_s = A_s * f_max_l * I_a / (0.5 * d)

        # TODO: Secondary side stringer contribution
        # Calculated in the same manner as primary side stringers
        M_sl = 0.0

        # Moment reacted by longerons
        M_l = 0.5 * M_ext - M_c - M_s - M_sl

        # Longeron area to resist this moment
        A_l = M_l * 0.5 * d / (f_max_l * self.i_xx / self.area)

        # Now the compression sector is evaluated
        # Effectiveness of cover in compression is based on Peery curved panel buckling
        # F_CCR = ( 9*(t_c/R)**(5/3) + 0.16*(t_c/L)**(1.3) + K_c*np.pi**2/(12*(1-nu_c**2)) ) * cover_e

        # ... needs more code ... #

        #

        return A_l

    def forced_cripping(self) -> None:
        """Not implmented yet. Will be the same as MinorFrames."""
        raise NotImplementedError(
            "This will follow similar procedures at the MinorFrame."
        )


@dataclass
class Bulkhead(Component):
    """Fuselage pressure bulkhead component.

    Pressure bulkhead are located at structural synthesis cuts by the user-
    determined input. Assumptions used in this approach are:
        1.  Construction is stiffened sheet design simply supported around
            the periphery.
        2.  Strip theory provides an adequate definition of maximum bending
            moment.
        3.  Stiffeners are of constant cross section basedon the maxium bending
            moment, at equal spacings, and oriented parallel to the shortest
            bulkhead dimension.
        4.  Web thickness base don maximum pressure is constant throughout the
            bulkhead surface.
        5.  Minor frame material is used for bulkhead construction.

    There are 3 assumed pressure loading types: Uniform, Triangular, and
    Trapezoidal. Uniform loading occurs from cabin or equipment compartment
    pressurization. Triangular occurs from fuel head. Trapezoidal loadings
    result from the combinations of fuel head and vent pressure. Fuel pressure
    results from the vehicle maneuver such that both positive and negative
    maneuvers are examined.
    """

    duf: float
    """Design Ultimate Factor.

    (2.0 for personnel environment, 1.5 for equipment)
    """

    K_r: float
    """Fatigue reduction factor (percent of Ftu)"""

    p_1: float
    """Triangular pressure."""

    p_2: float
    """Uniform pressure."""

    L: float
    """Height of pressurized surface."""

    t_w: float
    """Web field thickness.

    (Webs are assumed to be milled.)
    """

    t_l: float
    """Web land thickness.

    (Webs are assumed to be milled.)
    """

    d: float
    """Stiffener spacing."""

    t_s: float
    """Stiffener web thickness."""

    H: float
    """Stiffener cap width.

    (An I-beam cap geometry.)
    """

    def allowable_tensile_stress(self, K_r: float) -> float:
        """The design allowable tensile stress.

        Args:
            K_r: Fatigue reduction factor (percent of Ftu).

        Returns:
            Float of design allowable tensile stress, F_t.
        """
        return min(self.material.F_tu / self.duf, K_r * self.material.F_tu)

    def max_bending_moment(self) -> float:
        """Bending moment per unit width for the pressure loading."""
        x = -self.L * self.p_2 + self.L * np.sqrt(
            self.p_2**2 + self.p_2 * self.p_1 + (self.p_1**2) / 3
        ) * self.p_1 ** (-1)
        k = x / self.L
        M_max = (
            self.p_2 * self.L**2 * 0.5 * (k - k**2)
            + self.p_1 * self.L**2 * (k - k**3) / 6
        )

        return float(M_max)

    def stiffener_area(self, t_s: float, H: float) -> float:
        """Area of stiffener, including effective web.

        Args:
            t_s: stiffener thickness
            H: stiffener cap width

        Returns:
            Float of calculated area.
        """
        return 6 * t_s * H

    def stiffener_inertia(self, t_s: float, H: float) -> float:
        """Second moment of area of stiffener, including effective web.

        Second order terms of thickness are assumed to be negligible.

        Args:
            t_s: stiffener thickness
            H: stiffener cap width

        Returns:
            Float of second moment of area.
        """
        return 14 / 3 * t_s * H**3

    def web_thickness(self, d: float) -> Tuple[float, float]:
        """Evaluate web thickness.

        Web sizing based on combined bending and diaphragm action between
        stiffeners. Webs are assumed to be milled between supports. Thicknesses
        are derived by curve-fit approximation of theoretical plots (same
        analytical method as Cover diaphragm sizing).

        Args:
            d: stiffener spacing

        Returns:
            (t_w, t_l): Web field thickness, Web land thickness
        """
        t_w = (
            1.3769
            * d
            * (self.p_1 + self.p_2) ** 2.484
            * self.material.E**1.984
            / self.allowable_tensile_stress(self.K_r) ** 4.467
        )

        t_l = (
            1.646
            * d
            * (self.p_1 + self.p_2) ** 0.894
            * self.material.E**0.394
            / self.allowable_tensile_stress(self.K_r) ** 1.288
        )

        return (t_w, t_l)

    def stiffener_spacing(
        self,
        d_1: float = 2.0,
        d_2: float = 12.0,
        H_1: float = 1.0,
        H_2: float = 5.0,
        t_s1: float = 0.025,
    ) -> Tuple[float, float, float]:
        """Stiffener spacing optimization routine.

        Stiffener spacing search is initiated at minimum spacing and
        continued at fixed increments where three values of effective thickness
        are obtained. A three-point curve fit solution is used to determine
        the optimum spacing.

        Args:
            d_1: spacing lower bound
            d_2: spacing upper bound
            H_1: width lower bound
            H_2: width upper bound
            t_s1: min gauge thickness

        Returns:
            t_s: minimum stiffener thickness
            d: stiffener spacing for minimum thickness
            H: stiffener width for minimum thickness

        Raises:
            ValueError: Converged thickness does not result in spacing within bounds.
        """
        x = np.linspace(d_1, d_2, num=5)
        y = []

        # Flange crippling set to compressive yield stress
        # fcc = 0.312*np.sqrt(self.material.F_cy*self.material.E_c)*(4*self.t_s/H_i)**(3/4)
        B = self.material.F_cy / (
            0.312 * np.sqrt(self.material.F_cy * self.material.E_c)
        )

        for d_i in x:
            t_w, t_l = self.web_thickness(d_i)

            t_s = (
                3
                * self.max_bending_moment()
                * B ** (8 / 3)
                / (224 * self.material.F_cy)
            ) ** (1 / 3)
            H = 4 * t_s / B ** (4 / 3)
            t_bar = (
                t_w
                + (1.1 * H * (t_l - t_w) / d_i)
                + ((4 * H * t_s + H * (2 * t_s - t_l)) / d_i)
            )
            y.append(t_bar)
            print(f"{d_i:>4.1f}, {t_s:.4f}, {t_w:.4f}, {t_bar:.3f}")

        # calculate polynomial
        z = np.polyfit(x, y, 3)
        f = np.poly1d(z)

        # minimum equivalent thickness
        res = minimize_scalar(f, bounds=(d_1, d_2), method="bounded")
        # The minimum t_bar
        t_bar = res.fun
        d = res.x
        t_s = t_s1 if t_s <= t_s1 else t_s
        H = 4 * t_s / B ** (4 / 3)

        if H < H_1:
            raise ValueError(
                f"Stiffener width of {H:.1f} required for min weight (t_s={t_s:0.3f}), but is outside bounds [{H_1:.1f}, {H_2:.1}]"
            )
        elif H > H_2:
            raise ValueError(
                f"Stiffener width of {H:.1f} required for min weight (t_s={t_s:0.3f}), but is outside bounds [{H_1:.1f}, {H_2:.1}]"
            )

        return (t_s, d, H)


@dataclass
class Fuselage:
    """A Fuselage assembly.

    The fuselage structural components serve a wide range of functions. For the purposes of weight
    estimating and accounting, these structural components are categorized as either basic or
    secondary structure according to the definitions in MIL-STD-1374. The weight estimating
    approach is based on calculating weights at the line item level of the detail weight
    statement report form.

    The program estimates basic structure weight by sizing structural members to strength,
    stiffness, fatigue, and manufacturing requirements. These requirements are established
    through the analysis of design criteria, engineering data, and vehicle geometry.
    As such, there are a SIGNIFICANT number of inputs involved.

    The approach to sizing shell structure (cover, minor frames, longerons or stringers)
    is that of a multistation analysis. Bulkheads and major frames are sized to their
    individual load requirements. The weights of these basic structure elements are
    sensitive to factors such as geometry, type of construction, material properties,
    temperature, loads (and loading criteria), acoustic fatigue, local panel flutter,
    cutout size and location, stiffness requirements, and manufacturing limitations.

    Secondary structure component weight are estimated by rule-of-thumb and empirical methods.
    The weights of these items are sensitive to factors such as vehicle type and usage,
    design criteria, specific item function, and dimensional data.

    For multistation analysis: The external shell sectional geometry is represented as
    a family of shapes (rounded rectangles). External geometry is described at the nose,
    tail, and 8 intermediate stations. Shell structure is evaluate at a maximum of 19
    synthesis cuts, for which geometry is determined by interpolation between the
    described stations.
    """

    nose: Station
    """Geometry definition at the Nose Station."""

    tail: Station
    """Geometry definition at the Tail Station."""
