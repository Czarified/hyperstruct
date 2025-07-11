"""Test cases for the fuselage module."""

import pytest

from hyperstruct import Material
from hyperstruct.fuselage import Cover
from hyperstruct.fuselage import ForcedCrippling


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
