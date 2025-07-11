"""Test cases for the MajorFrame and Bulkhead classes."""

import numpy as np
import pandas as pd
import pytest
import pytest_check as pycheck

from hyperstruct import Material
from hyperstruct import Station
from hyperstruct.fuselage import Bulkhead
from hyperstruct.fuselage import MajorFrame


@pytest.fixture
def aluminum() -> Material:
    """Basic aluminum."""
    material = Material(
        rho=0.101,
        E=10.5e6,
        E_c=10.6e6,
        nu=0.33,
        F_tu=79.0e3,
        F_ty=42.1e3,
        F_cy=48.3e3,
        F_su=47.0e3,
        F_bru=142.0e3,
        F_bry=104.0e3,
        F_en=20.0e3,
        db_r=116,
    )
    return material


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


def test_frame(aluminum: Material, rounded: Station) -> None:
    """Test a major frame."""
    frame = MajorFrame(
        material=aluminum,
        fs_loc=rounded.number,
        loads=np.array(
            [
                # y, z, V, H, M
                [15.0, 56.5, -20876.0, 0.0, 0.0],
                [-15.0, 56.5, -20876.0, 0.0, 0.0],
            ]
        ),
        geom=rounded,
        fd=6.0,
    )
    # Verify we don't have the attribute before synthesis
    pycheck.is_false(hasattr(frame, "weight"))
    frame.synthesis()
    # Check against expected weight after synthesis
    pycheck.almost_equal(frame.weight, 60.8, rel=0.1)


def test_bulkhead(aluminum: Material) -> None:
    """Check a basic Bulkhead for functionality."""
    bulkhead = Bulkhead(material=aluminum, duf=2.0, K_r=0.6, p_1=3.2, p_2=7.1, L=30.0)
    # Verify we don't have the attribute before synthesis
    pycheck.is_false(hasattr(bulkhead, "weight"))
    bulkhead.synthesis()
    # Check against expected weight after synthesis
    pycheck.almost_equal(bulkhead.weight, 6, rel=0.1)


if __name__ == "__main__":  # pragma: no cover
    # This snippet just verifies the show method and get_coords.
    material = Material(
        rho=0.101,
        E=10.5e6,
        E_c=10.6e6,
        nu=0.33,
        F_tu=79.0e3,
        F_ty=42.1e3,
        F_cy=48.3e3,
        F_su=47.0e3,
        F_bru=142.0e3,
        F_bry=104.0e3,
        F_en=20.0e3,
        db_r=116,
    )

    obj = Station(
        orientation="FS",
        name="Rounded Rectangle",
        number=1728.0,
        width=60.0,
        depth=50.0,
        vertical_centroid=30.0,
        radius=24.61,
    )
    cir = Station(
        orientation="FS",
        name="Circle",
        number=1700.0,
        width=50.0,
        depth=50.0,
        vertical_centroid=30.0,
        radius=25.0,
    )

    frm = MajorFrame(
        material=material,
        fs_loc=cir.number,
        loads=np.array(
            [
                # y, z, V, H, M
                [15.0, 56.5, -20876.0, 0.0, 0.0],
                [-15.0, 56.5, -20876.0, 0.0, 0.0],
            ]
        ),
        geom=cir,
        fd=6.0,
    )
    print(f"Upper Panel = {frm.geom.upper_panel:.3f}")
    print(f" Side Panel = {frm.geom.side_panel:.3f}")
    print(f"Lower Panel = {frm.geom.lower_panel:.3f}")

    # Check coords
    # y, z = frm.geom.get_coords(np.radians(130), debug=True)
    # print(f"In quadrant 4: y={y:.1f}, z={z:.1f}\n")
    # coords = [(y, z)]

    # y, z = frm.geom.get_coords(np.radians(225), debug=True)
    # print(f"In quadrant 3: y={y:.1f}, z={z:.1f}")
    # coords = [(y, z)]

    num = 25

    frm.synthesis(num)
    print("Geometry Table")
    print(frm.cuts)

    sizing = pd.DataFrame(frm.results, columns=["w_j", "tcap", "t_w_str", "t_w_res"])
    print("Sizing Results Table")
    print(sizing)
    print(f"Frame weight = {frm.weight:.1f} [lbs]")

    _ = frm.show(show_coords=True, save=False)

    # Sweep over number of cuts to observe the prediction sensitivity
    # nums = np.linspace(15, 90, dtype=int, num=25)
    # weights = []
    # for n in nums:
    #     frm.synthesis(n)
    #     weights.append(copy(frm.weight))

    # ser = pd.Series(weights, index=nums, name="Frame Weight Sweep")
    # ser.index.name = "Cuts"
    # print(ser)

    # ser.plot(kind="bar", ylabel="weight")
    # plt.show()
