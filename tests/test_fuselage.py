"""Test cases for the fuselage module."""

from typing import Tuple

import numpy as np
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import ArrayLike

from hyperstruct import Material
from hyperstruct import Station
from hyperstruct.fuselage import Cover
from hyperstruct.fuselage import ForcedCrippling
from hyperstruct.fuselage import Fuselage
from hyperstruct.fuselage import MajorFrame


@pytest.fixture
def aluminum() -> Material:
    """Some basic aluminum."""
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
    return material


@pytest.fixture
def unmilled_cover(aluminum: Material) -> Cover:
    """Build a Cover component."""
    component = Cover(
        material=aluminum, milled=False, L=30, D=20, R=1, RC=25, V=420.0, I=69.0, Q=69.0
    )
    return component


@pytest.fixture
def diag_ten(aluminum: Material) -> ForcedCrippling:
    """Build a ForcedCrippling class."""
    component = ForcedCrippling(
        d=15,
        h=7.25,
        c=3.0,
        b=1.0,
        construction="stringer",
        frame_material=aluminum,
        cover_material=aluminum,
        long_material=aluminum,
        t_r=0.040,
    )
    return component


@pytest.fixture
def a_station() -> Tuple[Station]:
    """Placeholder Station."""
    return (
        Station(
            orientation="FS",
            name="Placeholder",
            number=200.0,
            width=40.0,
            depth=40.0,
            vertical_centroid=20.0,
            radius=15.5,
        ),
    )


@pytest.fixture
def b_station() -> Tuple[Station]:
    """Placeholder Station."""
    return (
        Station(
            orientation="FS",
            name="Placeholder",
            number=420.0,
            width=40.0,
            depth=40.0,
            vertical_centroid=20.0,
            radius=15.5,
        ),
    )


@pytest.fixture
def b_frame(aluminum: Material, b_station: Station) -> Tuple[MajorFrame]:
    """Placeholder MajorFrame."""
    return (
        MajorFrame(material=aluminum, fs_loc=420.0, loads=0.0, geom=b_station, fd=4.2),
    )


@pytest.fixture
def fuselage(a_station: Tuple[Station], b_frame: Tuple[MajorFrame]) -> Fuselage:
    """Build a Fuselage class."""
    fuse = Fuselage(stations=a_station, major_frames=b_frame)
    return fuse


@pytest.fixture
def fuse_loads() -> Tuple[ArrayLike, ArrayLike, ArrayLike, ArrayLike]:
    """Some basic beam loads."""
    w_fus = np.array(
        [
            [100.0, -245.0, 0.0],
            [435.2, -300.0, 0.0],
            [521.0, -450.0, 0.0],
        ]
    )
    w_fc = np.array(
        [
            [110.0, -200.0, 0.0],
            [400.5, -120.0, 0.0],
        ]
    )
    p_air = np.array(
        [
            [110.0, 200.0, 0.0],
            [400.5, 120.0, 0.0],
            [435.2, 300.0, 0.0],
        ]
    )
    p_ext = np.array(
        [
            [100.0, 245.0, 0.0],
            [521.0, 450.0, 0.0],
        ]
    )
    return (w_fus, w_fc, p_air, p_ext)


def test_unmilled_shear_and_net(unmilled_cover: Cover) -> None:
    """Test an unmilled cover."""
    t_c = unmilled_cover.field_thickness_block_shear()
    t_l = unmilled_cover.land_thickness_net_section()
    assert isinstance(t_l, float)
    assert isinstance(t_c, float)


def test_unmilled_pressure(unmilled_cover: Cover) -> None:
    """Test an unmilled cover."""
    t_l, t_c = unmilled_cover.thickness_pressure()
    assert isinstance(t_l, float)
    assert isinstance(t_c, float)


def test_unmilled_flutter(unmilled_cover: Cover) -> None:
    """Test an unmilled cover."""
    t_c = unmilled_cover.panel_flutter(mach=1.3, altitude=5000)
    assert isinstance(t_c, float)


def test_unmilled_acoustic(unmilled_cover: Cover) -> None:
    """Test an unmilled cover."""
    t_l, t_c = unmilled_cover.acoustic_fatigue()
    assert isinstance(t_l, float)
    assert isinstance(t_c, float)


def test_diagonal_tension(diag_ten: ForcedCrippling) -> None:
    """Initializes a ForcedCrippling class."""
    t_r, f_st, f_rg = diag_ten.forced_crippling(
        D=65.8, M=1.475e6, Z=32.9, sum_z_sq=120, t_c=0.032, RC=30.0, f_s=6000, f_scr=815
    )
    assert t_r <= 0.032
    assert f_st >= 15000
    assert f_rg >= 6500


def test_net_loads(fuselage: Fuselage, fuse_loads: Tuple[ArrayLike]) -> None:
    """Tests internal loads calculations on a Fuselage class."""
    w_fus, w_fc, p_air, p_ext = fuse_loads
    internal_loads = fuselage.net_loads(w_fus, w_fc, p_air, p_ext)
    assert internal_loads is not None


def test_vmt_diagram(fuselage: Fuselage, fuse_loads: Tuple[ArrayLike]) -> None:
    """Tests the diagram method."""
    w_fus, w_fc, p_air, p_ext = fuse_loads
    fig, axs = fuselage.vmt_diagram(w_fus, w_fc, p_air, p_ext)
    assert isinstance(fig, Figure)
    assert isinstance(axs[0], Axes)


if __name__ == "__main__":  # pragma: no cover
    # Plot the example fuse loads for visual check
    # Prototype the VMT diagram code
    import matplotlib.pyplot as plt

    w_fus = np.array(
        [
            [100.0, -245.0, 0.0],
            [435.2, -300.0, 0.0],
            [521.0, -450.0, 0.0],
        ]
    )
    w_fc = np.array(
        [
            [110.0, -200.0, 0.0],
            [400.5, -120.0, 0.0],
        ]
    )
    p_air = np.array(
        [
            [110.0, 200.0, 0.0],
            [400.5, 120.0, 0.0],
            [435.2, 300.0, 0.0],
        ]
    )
    p_ext = np.array(
        [
            [100.0, 245.0, 0.0],
            [521.0, 450.0, 0.0],
        ]
    )

    mat = Material(
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

    a = Station(
        orientation="FS",
        name="Placeholder",
        number=200.0,
        width=40.0,
        depth=40.0,
        vertical_centroid=20.0,
        radius=15.5,
    )
    b = Station(
        orientation="FS",
        name="Placeholder",
        number=200.0,
        width=40.0,
        depth=40.0,
        vertical_centroid=20.0,
        radius=15.5,
    )

    stations = (a,)
    frames = (MajorFrame(material=mat, fs_loc=420.0, loads=0.0, geom=b, fd=4.2),)

    fuse = Fuselage(stations=stations, major_frames=frames)
    loads = fuse.net_loads(w_fus, w_fc, p_air, p_ext)
    fig, (ax1, ax2) = fuse.vmt_diagram(w_fus, w_fc, p_air, p_ext)

    print("   FS     , P     , M_ext   , V     ,  M_int")
    print(loads)

    print("\n\n")
    x, v, m = fuse.lookup_loads(x=100, loads=loads)
    print(x)
    print(v)
    print(m)

    _ = ax1.plot(x, v, marker="^", color="k")
    _ = ax2.plot(x, m, marker="^", color="k")

    plt.show()
