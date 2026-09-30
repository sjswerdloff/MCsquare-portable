"""Water-cube CT and single-layer square-field plans for the #16 field-edge matrix.

Writes cube.mhd/cube.raw (water, 2 mm voxels, FE_N per side, default 150 = 300 mm) and one plan per (energy, field side),
named E<energy>_S<side mm>.txt, into the current directory. Prints the case names, one per line.
"""

import os
import sys

import numpy as np

N = int(os.environ.get("FE_N", "150"))  # voxels per side; 150 = a 300 mm cube
SP = 2.0
ISO = (N * SP / 2,) * 3
SPOT_SPACING = 5.0


def write_cube() -> None:
    """Write the water-cube CT (0 HU everywhere)."""
    with open("cube.mhd", "w") as f:
        f.write(
            f"ObjectType = Image\nNDims = 3\nDimSize = {N} {N} {N}\nElementSpacing = {SP:f} {SP:f} {SP:f}\n"
            "Offset = 0.000000 0.000000 0.000000\nElementType = MET_FLOAT\nElementByteOrderMSB = False\n"
            "ElementDataFile = cube.raw\n"
        )
    np.zeros((N, N, N), dtype=np.float32).tofile("cube.raw")


def write_plan(name: str, energy: float, side_mm: float) -> None:
    """Write a one-layer square field of equally weighted spots, range shifter OUT (as in the #16 test)."""
    half = side_mm / 2
    grid = np.arange(-half + SPOT_SPACING / 2, half, SPOT_SPACING)
    spots = [(x, y) for x in grid for y in grid]
    total = float(len(spots))
    lines = [
        "#TREATMENT-PLAN-DESCRIPTION", "#PlanName", name, "#NumberOfFractions", "1", "##FractionID", "1",
        "##NumberOfFields", "1", "###FieldsID", "1", "#TotalMetersetWeightOfAllFields", f"{total:f}", "",
        "#FIELD-DESCRIPTION", "###FieldID", "1", "###FinalCumulativeMeterSetWeight", f"{total:f}",
        "###GantryAngle", "0.000000", "###PatientSupportAngle", "0.000000",
        "###IsocenterPosition", f"{ISO[0]:f}\t {ISO[1]:f}\t {ISO[2]:f}",
        "###RangeShifterID", "RS_Block", "###RangeShifterType", "binary",
        "###NumberOfControlPoints", "1", "", "#SPOTS-DESCRIPTION",
        "####ControlPointIndex", "1", "####SpotTunnedID", "1", "####CumulativeMetersetWeight", "0.000000",
        "####Energy (MeV)", f"{energy:f}", "####RangeShifterSetting", "OUT",
        "####IsocenterToRangeShifterDistance", "300.000000", "####RangeShifterWaterEquivalentThickness", "74.100000",
        "####NbOfScannedSpots", str(len(spots)), "####X Y Weight",
    ]
    lines += [f"{x:f} {y:f} 1.000000" for x, y in spots]
    with open(f"{name}.txt", "w") as f:
        f.write("\n".join(lines) + "\n")


def main(energies: list[float], sides_mm: list[float]) -> None:
    write_cube()
    for e in energies:
        for s in sides_mm:
            name = f"E{e:.0f}_S{s:.0f}"
            write_plan(name, e, s)
            print(name)


if __name__ == "__main__":
    energies = [float(v) for v in sys.argv[1].split(",")]
    sides = [float(v) for v in sys.argv[2].split(",")]
    main(energies, sides)
