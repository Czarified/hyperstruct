"""Test cases for the fuselage module."""

import pytest

from hyperstruct import Material
from hyperstruct.fuselage import Cover


@pytest.fixture
def aluminum() -> Material:
    """Some basic aluminum."""
    material = Material(
        E=10.5e6,
        E_c=10.6e6,
        nu=0.33,
        F_tu=64,
        F_ty=42.1,
        F_cy=48,
        F_su=41,
        F_bru=104,
        F_bry=89,
        F_en=20,
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
