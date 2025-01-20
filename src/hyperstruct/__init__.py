"""Hyperstruct."""

from dataclasses import dataclass
from importlib.metadata import version

import numpy as np
from scipy.special import ellipeinc


__version__ = version("hyperstruct")


@dataclass
class Material:
    """The HyperStruct material model.

    The material object contains the basic material allowables.
    It also contains thermomechanical curves for environmental correction factors.
    """

    rho: float
    """Weight density (aka specific weight)"""

    E: float
    """Young's Modulus, tension."""

    E_c: float
    """Young's Modulus, compression."""

    nu: float
    """Poisson's ratio."""

    F_tu: float
    """Tensile ultimate strength."""

    F_ty: float
    """Tensile yield strength."""

    F_cy: float
    """Compressive yield strength."""

    F_su: float
    """Shear ultimate strength."""

    F_bru: float
    """Bearing ultimate strength."""

    F_bry: float
    """Bearing yield strength."""

    F_en: float
    """Endurance limit."""

    db_r: float
    """Random spectrum decibel level that produces an acoustic fatigue life of 10^9."""


@dataclass
class Component:
    """The basic foundation for sizing.

    An Aircraft is composed of different Assemblies, which contain various Components.
    Components have several base methods and attributes that enable sizing and
    weights estimation.
    """

    material: Material
    """material the cover is made of."""

    def synthesis(self) -> None:
        """The sizing method.

        The sizing method collects all sizing routines and executes them
        in the order of the `routines` list.
        """
        # This doesn't work. It's just a placeholder.
        pass


@dataclass
class Station:
    """The basic foundation for cross-sectional geometry.

    The station class provides functionality to define a cross-sectional shape,
    and location of that shape. These definitions are needed for multi-station
    analysis methods, such as the fuselage shell structure.

    Each Station class has an orientation system (FS, BL, WL), name, station
    number, and geometry definition.

    Differing from the SWEEP convention, shape geometry here is required
    to have a depth, width, and radius input. The correction factor approach
    seemed overly complex. If the radius provided is zero, a rectangle is
    generated. If the radius provided is greater than either the depth
    or width, an ellipse is assumed with major radius as the larger of the
    two and minor radius as the smaller of the two. The depth and width
    of the shape are required to align with the aircraft global coordinates
    of Z (WL. Vertical direction, from origin pointing UP) and Y (BL,
    Horizontal Plane, from origin pointing RIGHT when seated in the cockpit),
    respectively.

    Generally, sizing is accomplished for four shell sectors representing the
    upper, lower, and two sides. A 45-degree angle is used to define the
    limits of these sectors.
    """

    # This may not be necessary. Useful for Ben's reference.
    orientation: str
    """Fuselage Station, Butt Line, or Water Line."""

    name: str
    """A reference designation for this Station (e.g. Nose/Tail/MLG Trunnion)."""

    number: float
    """The station number. Measurement from the vehicle origin to the station plane."""

    width: float
    """Width of the shape."""

    depth: float
    """Depth of the shape."""

    vertical_centroid: float
    """Distance from z reference plane to the centroid."""

    radius: float
    """Corner radius of the shape."""

    @property
    def is_ellipse(self) -> bool:
        """Checks whether the shape is an ellipse."""
        if (self.radius > self.width / 2) or (self.radius > self.depth / 2):
            # An ellipse with major diameter equal to the larger dimension
            return True
        else:
            return False

    @property
    def doo(self) -> float:
        """Incremental depth from centroid to corner tangency."""
        if self.radius == 0:
            return 0.0
        elif self.is_ellipse:
            return 0.0
        else:
            return 0.5 * (self.depth - 2 * self.radius)

    @property
    def wo(self) -> float:
        """Incremental width from centroid to corner tangency."""
        if self.radius == 0:
            return 0.0
        elif self.is_ellipse:
            return 0.0
        else:
            return 0.5 * (self.width - 2 * self.radius)

    @property
    def upper_panel(self) -> float:
        """Peripheral length of the upper panel cover element."""
        if self.radius == 0:
            return self.width
        elif self.is_ellipse:
            start = 3 * np.pi / 4
            end = np.pi / 4
            arc = self.elliptical_arc_length(start, end)
            return arc
        else:
            return 2 * self.wo + np.pi / 2 * self.radius

    @property
    def lower_panel(self) -> float:
        """Peripheral length of the lower panel cover element."""
        if self.radius == 0:
            return self.width
        elif self.is_ellipse:
            start = -3 * np.pi / 4
            end = -np.pi / 4
            arc = self.elliptical_arc_length(start, end)
            return arc
        else:
            return 2 * self.wo + np.pi / 2 * self.radius

    @property
    def side_panel(self) -> float:
        """Peripheral length of the side panel cover elements."""
        if self.radius == 0:
            return self.depth
        elif self.is_ellipse:
            start = np.pi / 4
            end = -np.pi / 4
            arc = self.elliptical_arc_length(start, end)
            return arc
        else:
            return 2 * self.doo + np.pi / 2 * self.radius

    def elliptical_arc_length(self, theta_start: float, theta_end: float) -> float:
        """This method calculates the length of an elliptic curve.

        The semi-major axis is the distance from the centroid to the apex of
        the most extreme curve. The semi-minor axis is the distance from the
        centroid to the remaining apex.

        Note that the angles provided are measured in cartesian coordinates,
        counterclockwise, starting from the positive x-axis intersection of
        the curve.

        The semi-major axis is restricted to the larger of the Station depth
        or width, while the semi-minor axis is the opposite.

        Args:
            theta_start: Angle for start of curve, in radians (float).
            theta_end: Angle for end of curve, in radians (float).

        Returns:
            float: arc length
        """
        # semi-major radius
        a = max(self.width, self.depth)
        # semi-minor radius
        b = min(self.width, self.depth)
        # eccentricity
        m = 1 - (b / a) ** 2

        arc_length = a * (ellipeinc(theta_end, m) - ellipeinc(theta_start, m))

        return float(arc_length)
