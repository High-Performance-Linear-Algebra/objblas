#!/bin/bash
#
# Run GEMV roofline tests for different memory regimes.
# Output from each run is appended to results_gemv.csv (created by the program).
#

export OMP_NUM_THREADS=32
EXE=../../build/gemv_perf

# Machine constants for Intel i9-14900HX
BW_GBS=89.6
PEAK_GFLOPS=844.8

echo "Running GEMV roofline tests..."
echo "Executable: $EXE"
echo ""

run_test () {
    local M=$1
    local N=$2
    local NX=$3
    local NY=$4
    local REPS=$5
    local LABEL=$6

    echo "-----------------------------------------------------"
    echo "Running: $LABEL"
    echo "m=$M  n=$N  nx=$NX  ny=$NY  reps=$REPS"
    echo "-----------------------------------------------------"

    # Feed inputs to the executable using a here-document
    $EXE <<EOF
$BW_GBS
$PEAK_GFLOPS
$M
$N
$NX
$NY
$REPS
EOF

    echo ""
}

# ---------------------------------------------------------------------------
# L1-cache-resident test
# ---------------------------------------------------------------------------
# Very small tiles, hot in L1
run_test \
    256 256 \
    32 32 \
    2000 \
    "L1 test (small working set)"

# ---------------------------------------------------------------------------
# L2/L3-resident test
# ---------------------------------------------------------------------------
run_test \
    3000 3000 \
    256 256 \
    500 \
    "L2/L3 cache test (medium working set)"

# ---------------------------------------------------------------------------
# DRAM-bound test
# ---------------------------------------------------------------------------
run_test \
    20000 20000 \
    512 512 \
    75 \
    "DRAM test (large working set)"

echo "All tests finished."
