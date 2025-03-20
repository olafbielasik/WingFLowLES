#!/bin/bash
# runsimulation.sh - Script to run the full CFD LES simulation and visualisation pipeline

# Suspend execution in the event of an error
set -e

echo "========================================"
echo "I am starting a CFD simulation pipeline"
echo "========================================"

echo "Step 1: Generation of mesh and profile data (airfoil)"
python mesh_generation.py

echo "Step 2: Running the LES solver (results are saved every 500 steps)"
python solver_les.py

echo "Step 3: Generate a video animation from the saved frames"
python video_animation.py

echo "========================================"
echo "Simulation completed, the simulation.mp4 file has been generated."
echo "========================================"

