#!/bin/bash
# Test MCsquare arm64 build with sample data
# Runs 1M primaries and outputs to Outputs/

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Check binary exists
if [[ ! -x ./MCsquare_arm64 ]]; then
    echo "Error: MCsquare_arm64 not found. Run 'make MCsquare_arm64' first."
    exit 1
fi

# Create output directory
mkdir -p Outputs

# Create test config
cat > test_config.txt << 'EOF'
Num_Threads         0
Num_Primaries       1e6
CT_File             Sample_input_data/CT.mhd
HU_Density_Conversion_File   Scanners/default/HU_Density_Conversion.txt
HU_Material_Conversion_File  Scanners/default/HU_Material_Conversion.txt
BDL_Machine_Parameter_File   BDL/BDL_default_DN_RangeShifter.txt
BDL_Plan_File       Sample_input_data/PlanPencil.txt
Output_Directory    Outputs
Dose_MHD_Output     True
EOF

echo "Running MCsquare arm64 test with 1M primaries..."
./MCsquare_arm64 test_config.txt

# Verify output
if [[ -f Outputs/Dose.mhd && -f Outputs/Dose.raw ]]; then
    echo ""
    echo "Test PASSED - Output files generated:"
    ls -la Outputs/
else
    echo "Test FAILED - Output files not found"
    exit 1
fi
