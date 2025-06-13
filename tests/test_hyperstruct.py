"""Test cases for the general module."""

import numpy as np
import pytest
import pytest_check as pycheck

from hyperstruct import Material
from hyperstruct import Station
from hyperstruct.fuselage import MajorFrame


@pytest.fixture
def circle() -> Station:
    """Build a Station class that is circular."""
    obj = Station(
        orientation="FS",
        name="Circle",
        number=245.0,
        width=60.0,
        depth=60.0,
        vertical_centroid=60.0,
        radius=30.0,
    )
    return obj


@pytest.fixture
def oval_deep() -> Station:
    """Build a Station class that is oval."""
    obj = Station(
        orientation="FS",
        name="Deep",
        number=200.0,
        width=35.0,
        depth=60.0,
        vertical_centroid=60.0,
        radius=45.0,
    )
    return obj


@pytest.fixture
def oval_wide() -> Station:
    """Build a Station class that is oval."""
    obj = Station(
        orientation="FS",
        name="Wide",
        number=320.0,
        width=60.0,
        depth=35.0,
        vertical_centroid=60.0,
        radius=45.0,
    )
    return obj


@pytest.fixture
def rounded() -> Station:
    """Build a Station class that is rectangular."""
    obj = Station(
        orientation="FS",
        name="Rounded Rectangle",
        number=400.0,
        width=60.0,
        depth=60.0,
        vertical_centroid=60.0,
        radius=15.0,
    )
    return obj


def test_circle(circle: Station) -> None:
    """Test the Circular station."""
    upper = circle.upper_panel
    lower = circle.lower_panel
    side = circle.side_panel
    arc = circle.arc_length(0.0, np.pi / 4)
    x, y = circle.get_coords(np.pi / 2)

    pycheck.is_false(circle.is_ellipse)

    pycheck.almost_equal(upper, 47.1, abs=0.1)
    pycheck.almost_equal(lower, 47.1, abs=0.1)
    pycheck.almost_equal(side, 47.1, abs=0.1)

    pycheck.almost_equal(arc, 23.56, abs=0.1)

    pycheck.almost_equal(x, 30, abs=0.1)
    pycheck.almost_equal(y, 60, abs=0.1)


def test_oval_1(oval_deep: Station) -> None:
    """Test an Oval station."""
    upper = oval_deep.upper_panel
    lower = oval_deep.lower_panel
    side = oval_deep.side_panel
    arc = oval_deep.arc_length(0.0, np.pi / 4)
    x, y = oval_deep.get_coords(np.pi / 2)

    pycheck.is_true(oval_deep.is_ellipse)

    pycheck.almost_equal(upper, -31.7, abs=0.1)
    pycheck.almost_equal(lower, 31.7, abs=0.1)
    pycheck.almost_equal(side, -44.1, abs=0.1)

    pycheck.almost_equal(arc, 22.0, abs=0.1)

    pycheck.almost_equal(x, 51.4, abs=0.1)
    pycheck.almost_equal(y, 60, abs=0.1)


def test_oval_2(oval_wide: Station) -> None:
    """Test an Oval station."""
    upper = oval_wide.upper_panel
    lower = oval_wide.lower_panel
    side = oval_wide.side_panel
    arc = oval_wide.arc_length(0.0, np.pi / 4)
    x, y = oval_wide.get_coords(np.pi / 2)

    pycheck.is_true(oval_wide.is_ellipse)

    pycheck.almost_equal(upper, -31.7, abs=0.1)
    pycheck.almost_equal(lower, 31.7, abs=0.1)
    pycheck.almost_equal(side, -44.1, abs=0.1)

    pycheck.almost_equal(arc, 22.0, abs=0.1)

    pycheck.almost_equal(x, 30, abs=0.1)
    pycheck.almost_equal(y, 60, abs=0.1)


def test_rounded(rounded: Station) -> None:
    """Test a Rounded station."""
    x, y = rounded.get_coords(np.pi / 2)

    pycheck.is_false(rounded.is_ellipse)

    pycheck.almost_equal(x, 30, abs=0.1)
    pycheck.almost_equal(y, 60, abs=0.1)


if __name__ == "__main__":  # pragma: no cover
    # This snippet just verifies the show method and get_coords.
    material = Material(
        rho=0.1,
        E=10.5e6,
        E_c=10.6e6,
        nu=0.33,
        F_tu=64e3,
        F_ty=42.1e3,
        F_cy=48.3e3,
        F_su=41.0e3,
        F_bru=10.04e3,
        F_bry=89.0e3,
        F_en=20.0e3,
        db_r=116,
    )

    obj = Station(
        orientation="FS",
        name="Rounded Rectangle",
        number=400.0,
        width=60.0,
        depth=50.0,
        vertical_centroid=30.0,
        radius=24,
    )
    # angles = np.linspace(0, np.pi, num=10)
    # coords = []
    # for theta in angles:
    #     x, y = obj.get_coords(theta)
    #     coords.append((x, y))

    frm = MajorFrame(
        material=material,
        fs_loc=400.0,
        loads=np.array(
            [
                [15.0, 56.5, 200.0, 10.0, 0.0],
                [31.5, 30.0, 0.0, 0.0, 10.0],
                [-15.0, 56.5, 200.0, -10.0, 0.0],
                [-31.5, 30.0, 0.0, 0.0, -10.0],
            ]
        ),
        geom=obj,
        fd=6.4,
    )
    print(f"Upper Panel = {frm.geom.upper_panel:.3f}")
    print(f" Side Panel = {frm.geom.side_panel:.3f}")
    print(f"Lower Panel = {frm.geom.lower_panel:.3f}")

    # Check quadrants 2 and 3
    # y, z = frm.geom.get_coords(np.radians(300), debug=True)
    # print(f"In quadrant 2: y={y:.1f}, z={z:.1f}\n")
    # coords = [(y, z)]

    # y, z = frm.geom.get_coords(np.radians(225), debug=True)
    # print(f"In quadrant 3: y={y:.1f}, z={z:.1f}")
    # coords = [(y, z)]

    frm.synthesis(16)
    print(f"Frame weight = {frm.weight:.1f} [lbs]")

    zzf, inertias, cuts = frm.geometry_cuts(32)
    coords = [(row[5], row[6]) for row in cuts]
    # for row in coords:
    #     print(row)

    _ = frm.show(coords)
