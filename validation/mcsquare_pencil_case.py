"""MCsquare inputs for the stage-1 pencil beam (issue #32, validation/topas_design.md).

Writes into the current directory:
  water.mhd / water.raw  0 HU everywhere (read with Scanners/Water_Phantom), 1 mm voxels
  bdl_E<energy>.txt      two identical rows bracketing the energy, Mean = Nominal, single Gaussian,
                         spot sigma 3 mm, divergence 1e-6 rad (a floor), correlation 1e-3, energy spread 0
  plan_E<energy>.txt     one spot at (0, 0), range shifter absent

Geometry. MCsquare places the CT's first voxel corner at (0, 0, 0) and does not read the MHD Offset for
transport. Its internal axes are x (lateral), y (beam, gantry 0), z (lateral); the MHD DimSize is written in
that order, so DimSize = NX NY NZ with NY the depth. The isocentre is put on the shared corner of the four
central lateral voxels (x = NX/2, z = NZ/2 in mm). Measured on an asymmetric probe (60 x 120 x 80, 100 MeV):
at gantry 0 the beam travels toward -y and enters through the y = NY face, so the stage-1 case uses
--iso-y NY (the entrance surface) with --nozzle 1.

Usage:
  python mcsquare_pencil_case.py ENERGY [--nx 400 --ny 350 --nz 400] --iso-y MM [--nozzle MM]
"""

import argparse
import math

import numpy as np

SPOT_SIGMA_MM = 3.0
DIVERGENCE_RAD = 1e-6
CORRELATION = 1e-3
SM_DISTANCE_MM = 10000.0  # only enters for off-axis spots; nonzero by construction
ROW_HALF_GAP_MEV = 1.0


def write_ct(nx: int, ny: int, nz: int) -> None:
    """Write a water CT (0 HU) with 1 mm voxels, x fastest on disk."""
    with open("water.mhd", "w", newline="\n") as f:  # LF on every platform
        f.write(
            f"ObjectType = Image\nNDims = 3\nDimSize = {nx} {ny} {nz}\nElementSpacing = 1.000000 1.000000 1.000000\n"
            "Offset = 0.000000 0.000000 0.000000\nElementType = MET_FLOAT\nElementByteOrderMSB = False\n"
            "ElementDataFile = water.raw\n"
        )
    np.zeros((nz, ny, nx), dtype=np.float32).tofile("water.raw")


def write_bdl(energy: float, nozzle_mm: float) -> str:
    """Write the two-row BDL; returns its file name."""
    name = f"bdl_E{energy:.0f}.txt"
    s, d, c = SPOT_SIGMA_MM, DIVERGENCE_RAD, CORRELATION
    rows = []
    for e in (energy - ROW_HALF_GAP_MEV, energy + ROW_HALF_GAP_MEV):
        g = f"{s:f} \t {d:.9f} \t {c:f} \t {s:f} \t {d:.9f} \t {c:f}"
        rows.append(f"{e:.3f} \t {e:.3f} \t 0.000000 \t 1.0 \t 1.000000 \t {g} \t 0.000000 \t {g}")
    header = (
        "NominalEnergy \t MeanEnergy \t EnergySpread \t ProtonsMU \t Weight1 \t SpotSize1x \t Divergence1x \t "
        "Correlation1x \t SpotSize1y \t Divergence1y \t Correlation1y \t Weight2 \t SpotSize2x \t Divergence2x \t "
        "Correlation2x \t SpotSize2y \t Divergence2y \t Correlation2y"
    )
    text = [
        # Must match exactly: data_beam_model.c selects the UPenn reader on this line and otherwise falls
        # through, silently, to a legacy fscanf format that reads this file as garbage (every primary misses).
        "--UPenn beam model (double gaussian)--",
        "# two identical rows bracketing the energy; see validation/topas_design.md",
        "",
        "Nozzle exit to Isocenter distance",
        NOZZLE_FMT.format(nozzle_mm),  # NOZZLE_FMT: validate() checks this serialised value
        "",
        "SMX to Isocenter distance",
        f"{SM_DISTANCE_MM:.1f}",
        "",
        "SMY to Isocenter distance",
        f"{SM_DISTANCE_MM:.1f}",
        "",
        "Beam parameters",
        "2 energies",
        "",
        header,
        *rows,
    ]
    with open(name, "w", newline="\n") as f:  # LF on every platform
        f.write("\n".join(text) + "\n")
    return name


def write_plan(energy: float, iso: tuple[float, float, float]) -> str:
    """Write a one-spot plan at (0, 0); returns its file name."""
    name = f"plan_E{energy:.0f}.txt"
    lines = [
        "#TREATMENT-PLAN-DESCRIPTION",
        "#PlanName",
        f"pencil_E{energy:.0f}",
        "#NumberOfFractions",
        "1",
        "##FractionID",
        "1",
        "##NumberOfFields",
        "1",
        "###FieldsID",
        "1",
        "#TotalMetersetWeightOfAllFields",
        "1.000000",
        "",
        "#FIELD-DESCRIPTION",
        "###FieldID",
        "1",
        "###FinalCumulativeMeterSetWeight",
        "1.000000",
        "###GantryAngle",
        "0.000000",
        "###PatientSupportAngle",
        "0.000000",
        "###IsocenterPosition",
        "\t ".join(COORD_FMT.format(c) for c in iso),
        "###NumberOfControlPoints",
        "1",
        "",
        "#SPOTS-DESCRIPTION",
        "####ControlPointIndex",
        "1",
        "####SpotTunnedID",
        "1",
        "####CumulativeMetersetWeight",
        "0.000000",
        "####Energy (MeV)",
        f"{energy:f}",
        "####NbOfScannedSpots",
        "1",
        "####X Y Weight",
        "0.000000 0.000000 1.000000",
    ]
    with open(name, "w", newline="\n") as f:  # LF on every platform
        f.write("\n".join(lines) + "\n")
    return name


ENERGY_RANGE_MEV = (10, 300)
NOZZLE_FMT = "{:.1f}"  # as written into the BDL
COORD_FMT = "{:f}"  # as written into the plan


def serialised(nozzle: float, iso_y: float) -> tuple[float, float]:
    """The nozzle and isocentre depth exactly as the files will carry them (alden-ec2221c7, #40: 0.01 -> 0.0)."""
    return float(NOZZLE_FMT.format(nozzle)), float(COORD_FMT.format(iso_y))


def validate(energy: float, nx: int, ny: int, nz: int, iso_y: float, nozzle: float) -> None:
    """Refuse settings that would write a case MCsquare runs silently wrong, before anything is written.

    Raises:
        ValueError: naming the first violated condition.
    """
    for name, v in (("energy", energy), ("iso_y", iso_y), ("nozzle", nozzle)):
        if not math.isfinite(v):
            raise ValueError(f"{name} must be finite, got {v}")
    if energy != round(energy) or not ENERGY_RANGE_MEV[0] <= energy <= ENERGY_RANGE_MEV[1]:
        raise ValueError(f"energy must be a whole number of MeV in {ENERGY_RANGE_MEV} (file names use it), got {energy}")
    if min(nx, ny, nz) <= 0:
        raise ValueError(f"dimensions must be positive, got {nx} {ny} {nz}")
    if nx % 2 or nz % 2:
        raise ValueError("nx and nz must be even so the axis sits on the shared corner of the four central voxels")
    if not 0 <= iso_y <= ny:
        raise ValueError(f"iso_y must lie within the CT depth [0, {ny}] mm, got {iso_y}")
    if nozzle <= 0:
        raise ValueError(f"nozzle must be positive, got {nozzle}")
    if iso_y + nozzle <= ny:
        # The beam enters at y = NY travelling to -y; a source on or inside the CT starts in voxel index NY (one past
        # the grid) or mid-phantom, and the dose lands in the wrong place with exit 0.
        raise ValueError(f"source plane iso_y + nozzle = {iso_y + nozzle} mm must lie outside the CT (> {ny} mm)")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("energy", type=float)
    p.add_argument("--nx", type=int, default=400)
    p.add_argument("--ny", type=int, default=350, help="depth voxels (beam axis at gantry 0)")
    p.add_argument("--nz", type=int, default=400)
    p.add_argument("--iso-y", type=float, required=True, help="isocentre depth coordinate, mm")
    p.add_argument(
        "--nozzle",
        type=float,
        default=1.0,
        help="nozzle exit to isocentre, mm. NOT 0: with the isocentre on the entrance face a particle that "
        "starts there counts as inside the CT (Transport_to_CT's bounds are inclusive), sits at voxel "
        "index NY, one past the grid, and deposits everything in the first bin. 1 mm matches TOPAS's "
        "source 1 mm upstream; MCsquare's air polynomial then removes 0.47-0.77 keV (200-100 MeV).",
    )
    a = p.parse_args()
    try:
        for name, v in (("iso_y", a.iso_y), ("nozzle", a.nozzle)):
            if not math.isfinite(v):
                raise ValueError(f"{name} must be finite, got {v}")
        nozzle_s, iso_y_s = serialised(a.nozzle, a.iso_y)
        validate(a.energy, a.nx, a.ny, a.nz, iso_y_s, nozzle_s)  # the values the files will contain, not the inputs
    except ValueError as e:
        p.error(str(e))
    iso = (a.nx / 2, iso_y_s, a.nz / 2)
    write_ct(a.nx, a.ny, a.nz)
    print(write_bdl(a.energy, nozzle_s))
    print(write_plan(a.energy, iso))


if __name__ == "__main__":
    main()
