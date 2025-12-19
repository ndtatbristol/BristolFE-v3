# BristolFE v3
### Paul Wilcox
This repository contains a number of Matlab functions and example scripts for performing explicit time-marching Finite Element (FE) simulations for elastodynamic wave propagation. Simulations can be performed using the built-in solver that exploit's Matlab's native sparse-matrix and GPU processing capabilities or an alternative commercial solver such as [Pogo](www.pogo.software). Support for Abaqus may be added in the future.
## Installation
To use, clone (or download and unzip) the repository and either permanently add the `BristolFE-v3/code` folder to the Matlab path or include the line

`addpath(genpath('RELEVANT_PATH/BristolFE-v3/code'));`

at the start of any scripts that use the BristolFE functions.
## Getting started
The user interacts with BristolFE via the suite of documented Matlab functions in the folder `BristolFE-v3/code`. These provide tools for the creation, execution, and visualisation of models. They are intended to be called from user-written Matlab script files and a number of example scripts are provided in `BristolFE-v3/examples`. Most likely you will start with a copy of one of these and modify it according to your requirements.
## Summary of changes since v2
The v2 to v3 changes are designed to make the code more consistent and future-proof for expansion to 3D. The main changes are:
- Names of functions have been rationalised: functions exclusively for use in 2D models are prefixed `fn_2d_`; functions exclusively for use in 3D models are prefixed by `fn_2d_`; functions without these prefixes can be used in 2D or 3D models with the same arguments
- The entry point for all solvers is now `fn_FE_entry_point` (not `fn_BristolVE_v2`), with the solver being selected as an option (the default is `fe_options.solver = 'BristolFE'`).
- A cell array, `matls`, is now used to describe materials rather than a structure array. Each element in the model (except for interface elements) remains associated with a material via the `n_els x 1` vector, `mod.el_mat_i`, which contains integers referencing the corresponding cell in `matls`
- `el_types` is a new cell array of strings, one for each type of element used in the model
- `mod.el_typ_i` is no longer a cell array of strings describing element types, but instead `n_els x 1` vector of integer indices that reference the corresponding cell in `el_types`.
- `mod.el_abs_i` is still an`n_els x 1` vector of relative absorption of elements in the range 0 (no absorption) to 1 (maximum absorption) but is now mandatory in the `mod` structure, not optional

## Overview
The entry point function for solving a model is `res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options)`. 

When this function is called, a complete mesh must be specified (in mod), the materials used must be defined (in matls), the element types used must be defined (el_types) and one or more loading steps and the required outputs defined (in the cell array steps).

The requested results for the corresponding loading step are returned in the cell array res. Typical outputs are one or both of: 
    1. History outputs - complete time histories of the displacement (or pressure in fluids) at one or mode nodes, typically plotted as time-domain signals.
    2. Field output - snapshots of the complete wavefield (its local kinetic energy) at intervals in time, typically displayed as a movie and used as a visualisation tool.

Most of the code in the example scripts is concerned with preparing mod, matls, and steps before fn_FE_entry_point is called and then displaying the outputs.

EXAMPLES
========

In BristolFE-v3/examples you will find the following scripts which provide simple examples how to set up different features in models:
    1. fluid_example.m - simulate pressure waves in a fluid domain
    2. solid_example.m - simulate longitudinal and shear waves in a solid domain
    3. coupled_solid_fluid_example.m - simulate waves in a fluid domain coupled to a solid one, showing mode conversions at the interface
    4. absorbing_layer_example.m - same as 3 but this time with an absorbing layer on 3 sides of the domain to prevent reflections
    5. subdomain_example.m - simulate waves in a pristine domain, add a scatterer to a subdomain and combine results to obtain overall response
    6. subdomain_array_example.m - same as 5 but simulating FMC data from an array transducer
    7. solid_example_angled_excitation - same as 2 but with normal or shear forcing applied on angled edge of model to illustrate how to apply force at an angle


UPDATES SINCE PREVIOUS RELEASE
==============================

The main time-marching solver has been made more efficient. The previous issue of instability when using fluid-solid interface elements has been resolved by reformulating the solver to use velocity at the current time step rather than at the previous half-step. This is slightly less efficient per time step, but the lost in efficiency is more than compensated for by the unconditional stability, which enables time steps up to the CFL limit to be used without problems.


