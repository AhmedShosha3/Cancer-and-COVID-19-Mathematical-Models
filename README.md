# Cancer and COVID-19 Mathematical Models (MATLAB)

MATLAB simulations of a three-variable cancer model and an eight-compartment COVID-19 model, solved with integer-order (ODE45, Euler), fractional-order (Caputo-type Euler and a simplified PECE-style predictor-corrector) and variable-order methods.

![Cancer and COVID-19 Mathematical Models](Mathematical_Models.png)

## Contents

- `Cancer Model/`: tumor, healthy and immune cell model, with integer, fractional, variable-order and crossover simulations
- `Corona Model/`: eight-compartment COVID-19 model, with integer and fractional simulations and a real-data comparison
- `Results/` inside each model folder: selected output figures

## How to run

Open MATLAB, set the current folder to `Cancer Model/` or `Corona Model/`, and run a `Run_*.m` script. Each model folder is self-contained.

## Notes

- `Real_Data.m` compares the infected compartment with reported case data after scaling the data to the model's peak and shifting it in time. It is an alignment of shapes.
- Some scripts in `Cancer Model/` call the fractional solvers one time step at a time.
- Several scripts record videos or write files to the current folder.

## Team and contribution

Academic team graduation project.
My part: MATLAB implementation of the models and numerical methods, simulations and comparison plots.

## Author

Ahmed Ibrahim Mohamed Kamel · [GitHub](https://github.com/AhmedShosha3) · [LinkedIn](https://www.linkedin.com/in/ahmed-shosha-212192296/)
