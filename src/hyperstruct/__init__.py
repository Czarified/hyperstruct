"""Hyperstruct."""

from dataclasses import dataclass
from importlib.metadata import version
from typing import List
from typing import Optional
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import FancyBboxPatch
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
            arc = self.arc_length(start, end)
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
            arc = self.arc_length(start, end)
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
            arc = self.arc_length(start, end)
            return arc
        else:
            return 2 * self.doo + np.pi / 2 * self.radius

    def arc_length(self, theta_start: float, theta_end: float) -> float:
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

        arc_length = a * (ellipeinc(theta_end, m) - ellipeinc(theta_start, m)) / 2

        return float(arc_length)

    def curvature(self) -> float:
        """Calculates the nominal readius of curvature.

        Returns:
            float of the radius of curvature

        Raises:
            NotImplementedError: until developed
        """
        if True is False:
            return 0.0
        else:
            raise NotImplementedError

    def show(
        self, coords: Optional[List[Tuple[float, float]]] = None, display: bool = True
    ) -> Tuple:
        """Plot the station shape for a visual check.

        This method just uses matplotlib to draw the shape on a plot.
        It will select the appropriate shape (Artist) object, based on
        the Station properties, and put in a figure on it's own.
        """
        fig, ax = plt.subplots()
        lower_y = 0.0 if self.vertical_centroid >= 0 else self.vertical_centroid
        ax.set(
            xlim=(-1.1 * self.width, 1.1 * self.width),
            ylim=(lower_y, 1.1 * self.depth + self.vertical_centroid),
            aspect="equal",
        )

        if self.is_ellipse:
            obj = Ellipse(
                xy=(0, self.vertical_centroid),
                width=self.width,
                height=self.depth,
                fill=False,
                edgecolor="black",
                linewidth=1.5,
            )
        elif self.radius < self.width / 2:
            xy = (-self.width / 2, self.vertical_centroid - self.depth / 2)
            style = f"Round, pad=0.0, rounding_size={self.radius}"
            obj = FancyBboxPatch(
                xy=xy,
                width=self.width,
                height=self.depth,
                boxstyle=style,
                fill=False,
                edgecolor="black",
                linewidth=1.5,
            )
        else:
            obj = Circle(
                xy=(0, self.vertical_centroid),
                radius=self.radius,
                fill=False,
                edgecolor="black",
                linewidth=1.5,
            )
        ax.add_artist(obj)
        obj.set_clip_box(ax.bbox)

        # Plot the Corner Centers
        center_x = [self.wo, self.wo, -self.wo, -self.wo]
        center_y = [
            self.doo + self.vertical_centroid,
            -self.doo + self.vertical_centroid,
            -self.doo + self.vertical_centroid,
            self.doo + self.vertical_centroid,
        ]
        _ = ax.plot(
            center_x,
            center_y,
            marker="+",
            markerfacecolor="none",
            markeredgecolor="black",
            markersize=10,
            linestyle="none",
            alpha=0.8,
        )

        for x, y in zip(center_x, center_y, strict=True):
            _ = ax.add_artist(
                Circle(
                    (x, y),
                    radius=self.radius,
                    fill=False,
                    edgecolor="black",
                    linestyle=(0, (10, 14)),
                    linewidth=0.25,
                    alpha=0.7,
                )
            )

        if coords:
            for x, y in coords:
                ax.plot(
                    x,
                    y,
                    marker="o",
                    markerfacecolor="teal",
                    markeredgecolor="none",
                    markersize=4,
                )

        if display:
            plt.show()

        return (fig, ax)

    def _quadratic_sol(
        self, m: float, q: float, p: float
    ) -> Tuple[float, float, float, float]:
        """Breaking out the quadratic formula solution.

        Args:
            m: float of linear slope
            q: float of horizontal offset
            p: float of vertical offset

        Returns:
            tuple of the 2 intersections

        Raises:
            ValueError: if the curves don't intersect
        """
        # Plugging the linear equation into the circle equation and
        # solving for the dependent variable
        a = m**2 + 1
        b = 2 * (m * self.vertical_centroid - m * q - p)
        c = (
            q**2
            - self.radius**2
            + p**2
            - 2 * self.vertical_centroid * q
            + self.vertical_centroid**2
        )
        if b**2 - 4 * a * c < 0:
            # print(f"x={xx:.1f}, y={yy:.1f}")
            # print(f"theta={np.degrees(theta):.0f}deg, m={m:.1f}, does not intersect.")
            raise ValueError("Line does not intersect with shape!")

        x_1 = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
        x_2 = (-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
        y_1 = m * x_1 + self.vertical_centroid
        y_2 = m * x_2 + self.vertical_centroid

        return (x_1, x_2, y_1, y_2)

    def _corner_intersection(
        self,
        theta: float,
        phi_1: float,
        phi_2: float,
        phi_3: float,
        phi_4: float,
        phi_5: float,
        phi_6: float,
        phi_7: float,
        phi_8: float,
        debug: bool = False,
    ) -> Tuple[float, float, float, float]:
        """Just the corners.

        Args:
            theta: Angle of cut
            phi_1: Upper RH corner start
            phi_2: Upper RH corner end
            phi_3: Lower RH corner start
            phi_4: Lower RH corner end
            phi_5: Lower LH corner start
            phi_6: Lower LH corner end
            phi_7: Upper LH corner start
            phi_8: Upper LH corner end
            debug: print debug information

        Returns:
            tuple of the 2 intersections
        """
        # I think it's better to call all phi angles in order, but calling for
        # useless args will upset the CI gods. So we do something meaningless with them.
        __unused = phi_1 + phi_2 + phi_4 + phi_6 + phi_8
        __unused += __unused
        # Intersection is on corner
        # Set the equation for the line and the corner circle equal to each other
        # Slope of line, y-intercept is zero; y = mx + vertical_centroid
        if debug:
            print("Corner Intersection")
        m = np.cos(theta) / np.sin(theta)

        # The specific corner designates the circle offsets
        if theta < phi_3:
            # We're in the Upper RH
            if debug:
                print("Upper RH")
            p = self.wo
            q = self.vertical_centroid + self.doo
        elif theta < phi_5:
            # We're in the Lower RH
            if debug:
                print("Lower RH")
            p = self.wo
            q = self.vertical_centroid - self.doo
        elif theta < phi_7:
            # We're in the Lower LH
            if debug:
                print("Lower LH")
            p = -self.wo
            q = self.vertical_centroid - self.doo
        else:
            # We're in the Upper LH
            if debug:
                print("Upper LH")
            p = -self.wo
            q = self.vertical_centroid + self.doo

        x_1, x_2, y_1, y_2 = self._quadratic_sol(m, q, p)

        return (x_1, x_2, y_1, y_2)

    def _switch(
        self, x_1: float, x_2: float, y_1: float, y_2: float, debug: bool = False
    ) -> Tuple[float, float]:
        """Just the corners.

        Args:
            x_1: float
            x_2: float
            y_1: float
            y_2: float
            debug: bool

        Returns:
            tuple of the desired intersection point
        """
        if debug:
            print(f"    x_1 = {x_1:.2f}")
            print(f"    x_2 = {x_2:.2f}")
            print(f"    y_1 = {y_1:.2f}")
            print(f"    y_2 = {y_2:.2f}")
        # The intersection we are interested in is always the one with the largest
        # x-magnitude.
        x_i = np.array([x_1, x_2])
        y_i = np.array([y_1, y_2])
        x = x_i.flat[np.abs(x_i).argmax()]
        y = y_i[0] if x == x_i[0] else y_i[1]

        return (x, y)

    def _rect_corner(self, theta: float, debug: bool = False) -> Tuple[float, float]:
        """This method breaks out all the rounded rectangle stuff corner stuff from get_coords.

        Args:
            theta:  polar angle to calculate the intersection with. This must be in radians.
                    Note that the angle is measured clockwise from the vertical axis (12 o'clock)
            debug:  bool to print debugging information

        Returns:
            Tuple of coordinates as floats, for example (x, y)
        """
        # Rounded rectangle
        # Start with the inscribed circle xy
        r = self.depth / 2
        xx = r * np.sin(theta)
        yy = r * np.cos(theta) + self.vertical_centroid

        if debug:
            print(f"xx={xx:.1f}, yy={yy:.1f}")

        # Corner points to polar angles
        d = self.depth / 2
        w = self.width / 2

        phi_1 = np.arctan(self.wo / d)
        phi_2 = np.arctan(w / self.doo)
        phi_3 = np.arctan(self.doo / w) + np.pi / 2
        phi_4 = np.pi - np.arctan(self.wo / d)
        phi_5 = np.arctan(self.wo / d) + np.pi
        phi_6 = 3 * np.pi / 2 - np.arctan(self.doo / w)
        phi_7 = 3 * np.pi / 2 + np.arctan(self.doo / w)
        phi_8 = 3 * np.pi / 2 + np.arctan(d / self.wo)

        # Stupid complex switch case
        if (theta < phi_1) or (theta > phi_8):
            # Intersection is on upper horizontal portion
            x = r * np.tan(theta)
            y = r + self.vertical_centroid
            if debug:
                print(f"Upper Horizontal Intersection; x={x:.1f}, y={y:.1f}")
        elif (theta > phi_4) and (theta < phi_5):
            # Intersection is on lower horizontal portion
            x = -r * np.tan(theta)
            y = -r + self.vertical_centroid
            if debug:
                print(f"Lower Horizontal Intersection; x={x:.1f}, y={y:.1f}")
                print(f"theta = {np.degrees(theta):.2f}deg")
                print(f"phi_4 = {np.degrees(phi_4):.0f}deg")
                print(f"phi_5 = {np.degrees(phi_5):.0f}deg")
        elif (theta > phi_2) and (theta < phi_3):
            # Intersection is on RH vertical portion
            # Use a triangle instead of a circle
            x = w
            y = w * np.cos(theta) + self.vertical_centroid
            if debug:
                print(f"RH Vertical Intersection; x={x:.1f}, y={y:.1f}")
        elif (theta > phi_6) and (theta < phi_7):
            # Intersection is on LH vertical portion
            # Use a triangle instead of a circle
            x = -w
            y = w * np.cos(theta) + self.vertical_centroid
            if debug:
                print(f"LH Vertical Intersection; x={x:.1f}, y={y:.1f}")
        else:
            x_1, x_2, y_1, y_2 = self._corner_intersection(
                theta, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6, phi_7, phi_8, debug
            )
            x, y = self._switch(x_1, x_2, y_1, y_2, debug)

        return (float(x), float(y))

    def get_coords(self, theta: float, debug: bool = False) -> Tuple[float, float]:
        """Get the planar coordinates of the shape intersection with a straight line at angle theta.

        Coordinate extraction varies from the simplistic (Circle), to the very complext case of the
        rounded rectangle. The method used to evaluate the coordinates depends on the shape family
        (circle, ellipse, or rounded rectangle).

        Args:
            theta:  polar angle to calculate the intersection with. This must be in radians.
                    Note that the angle is measured clockwise from the vertical axis (12 o'clock)
            debug:  bool to print debugging information

        Returns:
            Tuple of coordinates as floats, for example (x, y)
        """
        if self.is_ellipse:
            # Calculate the eccentricity
            a = self.width / 2
            b = self.depth / 2
            major = max(a, b)
            minor = min(a, b)
            e = np.sqrt(1 - (minor / major) ** 2)
            alpha = np.pi / 2 - theta
            # Vector length to curve intersection point
            r = b / np.sqrt(1 - (e * np.cos(alpha)) ** 2)
            # Use trig to get the xy instead of polar
            x = r * np.sin(theta)
            y = r * np.cos(theta) + self.vertical_centroid
        elif self.radius < self.width / 2:
            x, y = self._rect_corner(theta, debug)
        else:
            # It's a circle
            r = self.radius
            # Use trig to get the xy instead of polar
            x = r * np.sin(theta)
            y = r * np.cos(theta) + self.vertical_centroid

        return (float(x), float(y))
