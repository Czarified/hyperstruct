"""Test cases for the fuselage module."""

import pytest

from hyperstruct import Material
from hyperstruct.fuselage import Bulkhead
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
def bulkhead(aluminum: Material) -> Bulkhead:
    """Build a Bulkhead component."""
    component = Bulkhead(
        material=aluminum,
        p_1=30.0,
        p_2=6.0,
        L=34.0,
        K_r=0.53,
        duf=1.5,
        t_w=0.027,
        t_l=0.080,
        d=2.0,
        t_s=0.090,
        H=1.5,
    )
    return component


@pytest.fixture
def diag_ten(aluminum: Material) -> ForcedCrippling:
    """Build a ForcedCrippling class."""
    component = ForcedCrippling(
        material=aluminum, c=4.0, b=2.0, construction="stringer", t_r=0.050
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


def test_bulkhead_stiffener_spacing(bulkhead: Bulkhead) -> None:
    """Test a pressure bulkhead."""
    print(bulkhead.allowable_tensile_stress(0.53))
    print(bulkhead.web_thickness(2.0))
    t_s, d, H = bulkhead.stiffener_spacing()
    assert t_s >= 0.025
    assert t_s < 0.500


def test_diagonal_tension(diag_ten: ForcedCrippling) -> None:
    """Initializes a ForcedCrippling class."""
    assert diag_ten.c >= 4.0
