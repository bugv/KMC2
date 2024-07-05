# Vacancy Mediated Substitutional Diffusion in Random Alloys

This library provides the diffusion of a single vacancy in an alloy using the kinetic Monte Carlo algorithm. 

Written by Zoë Evans (zoe.evans@epfl.ch) in a Spring 2024 Semester Project at the MADES EPFL lab.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)

## Installation

This code was developed in python 3.10.

It requires the following libraries:

- matplotlib
- numpy
- pandas
- pymatgen
- pytest


## Usage

### Run a single input file

The kinetic Monte Carlo diffusion of a single vacancy can be run for a structure that is specifically defined or a random alloy created by providing the desired composition of the alloy.
To run the code the KMC, the `main.py` file should be executes with either the `-exact` option or the `-binary` option.


```
python3 main.py -binary Input_binary.json
python3 main.py -exact Input.json
```

An example file for a specifically defined structure is given in `Input_Al_fcc.json` and an example of an input file to create a random alloy is given in `Input_binary.json`.

This will produce the following files:

 - initialized_struct.json :  a file containing the structure after the initialization.
 - POSCAR : A poscar file to visualize the structure
 - results.json : a file with the results after the KMC

Once the KMC has been done, it is possible to run `analysis.py` to calculate the diffusion coefficients.


### Run KMC for a range of concentrations

The kinetic Monte Carlo diffusion of a single vacancy in a random alloy at a range of concentrations can be done using the `range_binaries.py` script.
The first section of the script should be modified to choose the correct parameters.

This will produce the following files :

- an initialize structure json file for each concentration
- a results.dat file containing the L_ij parameters for each concentration
- a plot of the L_ij values as a function of the concentration