"""Test cases for the general module."""

import numpy as np
import pytest
import pytest_check as pycheck

from hyperstruct import Station


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

    pycheck.is_false(circle.is_ellipse)

    pycheck.almost_equal(upper, 47.1, abs=0.1)
    pycheck.almost_equal(lower, 47.1, abs=0.1)
    pycheck.almost_equal(side, 47.1, abs=0.1)

    pycheck.almost_equal(arc, 23.56, abs=0.1)


def test_oval_1(oval_deep: Station) -> None:
    """Test an Oval station."""
    upper = oval_deep.upper_panel
    lower = oval_deep.lower_panel
    side = oval_deep.side_panel
    arc = oval_deep.arc_length(0.0, np.pi / 4)

    pycheck.is_true(oval_deep.is_ellipse)

    pycheck.almost_equal(upper, -31.7, abs=0.1)
    pycheck.almost_equal(lower, 31.7, abs=0.1)
    pycheck.almost_equal(side, -44.1, abs=0.1)

    pycheck.almost_equal(arc, 22.0, abs=0.1)


def test_oval_2(oval_wide: Station) -> None:
    """Test an Oval station."""
    upper = oval_wide.upper_panel
    lower = oval_wide.lower_panel
    side = oval_wide.side_panel
    arc = oval_wide.arc_length(0.0, np.pi / 4)

    pycheck.is_true(oval_wide.is_ellipse)

    pycheck.almost_equal(upper, -31.7, abs=0.1)
    pycheck.almost_equal(lower, 31.7, abs=0.1)
    pycheck.almost_equal(side, -44.1, abs=0.1)

    pycheck.almost_equal(arc, 22.0, abs=0.1)


def test_rounded(rounded: Station) -> None:
    """Test a Rounded station."""
    upper = rounded.upper_panel
    lower = rounded.lower_panel
    side = rounded.side_panel
    arc = rounded.arc_length(0.0, np.pi / 4)

    pycheck.is_false(rounded.is_ellipse)

    # pycheck.almost_equal(upper, -31.7, abs=0.1)
    # pycheck.almost_equal(lower, 31.7, abs=0.1)
    # pycheck.almost_equal(side, -44.1, abs=0.1)

    # pycheck.almost_equal(arc, 22.0, abs=0.1)


if __name__ == "__main__":  # pragma: no cover
    # This is sample code from the Matplotlib Ellipse example.
    # Playground stuff to learn how to plot these shapes.
    # import matplotlib.pyplot as plt
    # import numpy as np
    # from matplotlib.patches import Ellipse

    # # Fixing random state for reproducibility
    # np.random.seed(19680801)

    # NUM = 10

    # ells = [
    #     Ellipse(
    #         xy=np.random.rand(2) * 10,
    #         width=np.random.rand(),
    #         height=np.random.rand(),
    #         angle=np.random.rand() * 360,
    #         fill=False,
    #         edgecolor=np.random.rand(3),
    #     )
    #     for i in range(NUM)
    # ]

    # fig, ax = plt.subplots()
    # ax.set(xlim=(0, 10), ylim=(0, 10), aspect="equal")

    # for e in ells:
    #     ax.add_artist(e)
    #     e.set_clip_box(ax.bbox)

    # plt.show()
    obj = Station(
        orientation="FS",
        name="Rounded Rectangle",
        number=400.0,
        width=60.0,
        depth=60.0,
        vertical_centroid=60.0,
        radius=20,
    )
    angles = np.radians([-15, 5, 15, 35, 50, 65, 80])
    coords = []
    for theta in angles:
        print(f"Angle = {theta:.2f}")
        x, y = obj.get_coords(theta)
        coords.append((x, y))

    obj.show(coords)
