"""Test cases for the general module."""

import pytest
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import FancyBboxPatch

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


if __name__ == "__main__":  # pragma: no cover
    # This is sample code from the Matplotlib Ellipse example.
    # Playground stuff to learn how to plot these shapes.
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Ellipse

    # Fixing random state for reproducibility
    np.random.seed(19680801)

    NUM = 10

    ells = [
        Ellipse(
            xy=np.random.rand(2) * 10,
            width=np.random.rand(),
            height=np.random.rand(),
            angle=np.random.rand() * 360,
            fill=False,
            edgecolor=np.random.rand(3),
        )
        for i in range(NUM)
    ]

    fig, ax = plt.subplots()
    ax.set(xlim=(0, 10), ylim=(0, 10), aspect="equal")

    for e in ells:
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        # e.set_alpha(np.random.rand())
        # e.set_facecolor(np.random.rand(3))

    plt.show()
