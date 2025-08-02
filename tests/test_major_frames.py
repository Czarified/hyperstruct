"""Test cases for the MajorFrame and Bulkhead classes."""

import numpy as np
import pandas as pd
import pytest
import pytest_check as pycheck

from copy import copy
import matplotlib.pyplot as plt

from hyperstruct import Material, composite_cg
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

    loads_list = [
        np.array( # FS420
            [
                # y, z, V, H, M
                [ 15.5, 26, -2500.0, 0.0, 32441.0],
                [-15.5, 26, -2500.0, 0.0,-32441.0],
                [ 5.0, 10,  1750.0, 0.0, 0.0],
                [-5.0, 10,  1750.0, 0.0, 0.0],
            ]
        ),
        np.array( # FS517
            [
                # y, z, V, H, M
                [ 15.5, 56.5, 10472.0, 0.0, 0.0],
                [-15.5, 56.5, 10472.0, 0.0, 0.0],
            ]
        ),
        np.array( # FS720
            [
                # y, z, V, H, M
                [ 18.0, 56.5, -20876.0, 0.0, 41.2e4],
                [-18.0, 56.5, -20876.0, 0.0, 41.2e4],
            ]
        ),
        np.array( # FS872
            [
                # y, z, V, H, M
                [ 2.0, 47.0, 0.0, 0.0, 1776.9],
                [ 2.0, 45.5, 0.0, 3421.0, 0.0],
            ]
        ),
    ]
    data ={
        420.00:((30, 32, 10), loads_list[0]),
        517.50:((60, 50, 20), loads_list[1]),
        720.00:((60, 50, 24), loads_list[2]),
        872.00:((30, 30, 14), loads_list[3])
    }

    stations = []
    frames = []
    for fs, (dims, load) in data.items():
        w,d,r = dims
        obj = Station(
                orientation="FS",
                name="Rounded Rectangle",
                number=fs,
                width=w,
                depth=d,
                vertical_centroid=30.0,
                radius=r,
            )
        stations.append(copy(obj))

        frames.append(
            MajorFrame(
                material=material,
                fs_loc=obj.number,
                loads=load,
                geom=obj,
                fd=6.0,
            )
        )


    num = 25
    weights = []
    for frm in frames:
        frm.synthesis(num)
        print(f"Frame FS{frm.fs_loc:.0f} weight = {frm.weight:.1f}[lbs]")
        weights.append((frm.weight, frm.fs_loc))
        _fig, ax = frm.show(show_coords=True, save=False)

    total_weight, cg = composite_cg(weights)
    print(f"Total Frame Weight = {total_weight:.2f}[lbs] 😬")
    print(f"CG of Frames = {cg:.2f}[in] ⚖")
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
