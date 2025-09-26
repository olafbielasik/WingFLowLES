# Wing Flow - Large eddy simulation model

This project utilizes a Large Eddy Simulation (LES) combined with the Immersed Boundary (IB) method to model turbulent flow around a Clark Y airfoil at a -20-degree angle of attack. The simulation analyzes flow field dynamics, including velocity components (u, v), pressure (p), and vorticity (ω), while producing visualizations such as vorticity plots and video animations. Notably, it highlights the lift force generation resulting from the pressure and velocity distribution around the airfoil.

<h2 align="center">⬇️ Watch on YouTube ⬇️</h2>

<p align="center">
  <a href="https://www.youtube.com/watch?v=CqgccimCQGE">
    <img src="https://img.youtube.com/vi/CqgccimCQGE/maxresdefault.jpg" alt="Turbulent Flow Simulation" width="800"/>
  </a>
</p>

## Project Overview

The simulation models turbulent flow over a Clark Y airfoil using the LES approach with the Smagorinsky subgrid-scale model. A non-uniform grid with enhanced resolution near the airfoil is employed to resolve fine-scale turbulence. Key features include:
- **Airfoil Geometry**: Clark Y airfoil scaled to a chord length of 1.5 units, rotated by -20 degrees
- **Domain**: Rectangular domain with x ∈ [-0.5, 3.0] and y ∈ [-1.0, 1.0], discretized into a 600x400 grid
- **Simulation Parameters**: Reynolds number (Re) = 2000, Smagorinsky constant (Cs) = 0.1, time step (dt) = 0.0002, total steps = 60,000
- **Outputs**: Velocity fields (u, v), pressure (p), vorticity (ω), saved every 500 steps, and a final video animation (`simulation.mp4`)

This simulation primarily relies on CPU computations for solving the filtered Navier-Stokes equations
The solver performs approximately 9.67 trillion operations over the entire simulation
For an estimate of computational operations, see `computation_estimate.txt`
For a mathematical description see `mathematical_descritpion.txt`
For the simulation video, the duration of 5 seconds (real-time) was extended to 12 seconds, with interpolation to 120 FPS using FlowFrames and upscaling to 4K resolution using Video2X

<p align="center">
  <img src="data/simulation.gif" alt="WingFlowLES Simulation">
</p>

## Requirements

To run the simulation, set up the environment using the provided `environment.yml` file. Then execute the full pipeline using the provided Bash script in `run_simulation.sh`

## License

This project is licensed under the MIT License - see the LICENSE.txt file for details. 



