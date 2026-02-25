# Validation Tools

Developer tools for comparing MCsquare dose distributions across builds, platforms, and code changes. **Not part of the MCsquare deliverable.**

## Quick Start

```bash
pip install -r validation/requirements.txt

# Compare two dose outputs (text summary)
python validation/compare_doses.py reference/Dose.mhd evaluation/Dose.mhd

# With plots
python validation/compare_doses.py reference/Dose.mhd evaluation/Dose.mhd --output report/

# Custom gamma criteria (2%/2mm)
python validation/compare_doses.py ref/Dose.mhd eval/Dose.mhd --criteria 2 2

# Skip gamma (no pymedphys dependency needed)
python validation/compare_doses.py ref/Dose.mhd eval/Dose.mhd --no-gamma
```

## What It Compares

- **Dose difference statistics**: max, mean, RMS differences (absolute and % of max)
- **Depth-dose curves (PDD)**: central axis comparison
- **Lateral profiles**: crossline profiles at specified depths
- **Gamma analysis**: using pymedphys (default 3%/3mm, 20% lower cutoff)

## Typical Workflow

```bash
# Run reference build (e.g., original Intel/MKL)
./MCsquare_linux config.txt
mv Outputs reference_output

# Run portable build
./MCsquare_portable config.txt
mv Outputs portable_output

# Compare
python validation/compare_doses.py \
    reference_output/Dose.mhd \
    portable_output/Dose.mhd \
    --output comparison_report/
```

## Pass Criteria

Gamma pass rate >= 95% at the specified criteria (default 3%/3mm).

Exit code 0 = PASS, exit code 1 = FAIL.
