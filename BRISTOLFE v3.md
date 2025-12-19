# BristolFE v3
### Paul Wilcox
This repository contains a number of Matlab functions and example scripts for performing explicit time-marching Finite Element (FE) simulations for elastodynamic wave propagation. Simulations can be performed using the built-in solver that exploit's Matlab's native sparse-matrix and GPU processing capabilities or an alternative commercial solver such as [Pogo](www.pogo.software). Support for Abaqus may be added in the future.
## Installation
To use, clone (or download and unzip) the repository and either permanently add the `BristolFE-v3/code` folder to the Matlab path or include the line

`addpath(genpath('RELEVANT_PATH/BristolFE-v3/code'));`

at the start of any scripts that use the BristolFE functions.
## Getting started
The user interacts with BristolFE via the suite of documented Matlab functions in the folder `BristolFE-v3/code`. These provide tools for the creation, execution, and visualisation of models. They are intended to be called from user-written Matlab script files and a number of example scripts are provided in `BristolFE-v3/examples`. Most likely you will start with a copy of one of these and modify it according to your requirements.
## Examples
In BristolFE-v3/examples you will find the following scripts which provide simple examples how to set up different features in models:
- ex_2d_basic.m - script demonstrates a model with the minimum complexity
- ex_2d_advanced.m - script shows how to use more advanced features (e.g. fluid-solid coupling; absorbing layers; irregular inclusions; zero-volume cracks)
- ex_2d_subdomain.m - script demonstrates how to use 2D subdomains for more efficient modelling of localised scatterers
- ex_2d_pogo.m - script runs a simple simulation twice, once with the built-in solved and once with Pogo and then overlays the results (requires Pogo to be installed)
- ex_3d_pogo.m - script demonstrates 3D model creation tools and uses the Pogo solver (requires Pogo to be installed)
- ex_3d_subdomain.m - script demonstrates how to use 3D subdomains for more efficient modelling of localised scatterers (requires Pogo to be installed)
## Summary of changes since v2
The v2 to v3 changes are designed to make the code more consistent and future-proof for expansion to 3D. The main changes are:
- The built-in solver has been made more efficient and the previous issue of instability when using fluid-solid interface elements has been resolved by reformulating the solver to use velocity at the current time step rather than at the previous half-step
- Names of functions have been rationalised: functions exclusively for use in 2D models are prefixed `fn_2d_`; functions exclusively for use in 3D models are prefixed by `fn_2d_`; functions without these prefixes can be used in 2D or 3D models with the same arguments
- The entry point for all solvers is now `fn_FE_entry_point` (not `fn_BristolVE_v2`), with the solver being selected as an option (the default is `fe_options.solver = 'BristolFE'`).
- A cell array, `matls`, is now used to describe materials rather than a structure array. Each element in the model (except for interface elements) remains associated with a material via the `n_els x 1` vector, `mod.el_mat_i`, which contains integers referencing the corresponding cell in `matls`
- `el_types` is a new cell array of strings, one for each type of element used in the model
- `mod.el_typ_i` is no longer a cell array of strings describing element types, but instead `n_els x 1` vector of integer indices that reference the corresponding cell in `el_types`.
- `mod.el_abs_i` is still an`n_els x 1` vector of relative absorption of elements in the range 0 (no absorption) to 1 (maximum absorption) but is now mandatory in the `mod` structure, not optional

## Overview
The entry point function for solving a model is `res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options)`. 

When this function is called, a complete mesh must have been specified (in the structure `mod`), the materials used must have been defined (in the cell array `matls`), the element types used must have been defined (in the cell array `el_types`) and one or more loading steps and the required outputs must have been defined (in the cell array `steps`).

### Model description (`mod`)
This describes the model geometry and must contain the following fields:
- `mod.nds' - an `n_nds x n_dims` matrix of coordinates of all nodes in the model (`n_nds` is number of nodes; `n_dims` is number of dimensions, i.e. 2 or 3)
- `mod.els` - an `n_els x max_nds_per_el` matrix of the nodes associated with each element in the model (`n_els` is number of elements; `max_nds_per_el` is the maximum number of nodes used by an element in the model, which is typically 3 for a 2D model of triangular elements but can be as many as 8 for a 3D model with hexahedral elements)
- `mod.el_mat_i` - an `n_els x 1` vector of indices that describe the material associated with each element
- `mod.el_typ_i` - an `n_els x 1` vector of indices that describe the element type associated with each element
- `mod.el_abs_i` - an `n_els x 1` vector of numbers in the range 0 to 1 that describe the relative damping of each element (damping can be specified directly as a material property, but the relative damping controls how the elements damping and stiffness behaviour is modified when it is part of an absorbing layer)

*Note that node and element numbers are implicitly defined by the associated row number of the relevant matrix, where the first row represents node or element 1 (not 0).*
 
 ### Material descriptions (`matls`)
 The materials to be used in a model are defined in a `matls` cell array, with the cell index being the identifier references by `mod.el_mat_i`. The requried fields are:
 - `matls{i}.name` - a string giving the name
 - `matls{i}.rho` - the density
 - `matls{i}.D` - either:
   - solids - a `6 x 6` stiffness matrix (in a 2D model the necessary reduction to a `3 x 3` stiffness matrix takes place when plane stress or plane straing elements are formed)
   - fluids - a `1 x 1` matrix containing the bulk modulus
 - `matls{i}.col` - an `1 x 3` vector of RBG values
Typically functions such as the following deal can convert engineering data into the necessary values for some common cases:
- `matls{i} = fn_matl_isotropic_solid_defined_by_velocities(name, longitudinal_velocity, shear_velocity, density)`
- `matls{i} = fn_matl_fluid_defined_by_velocity(name, velocity, density)`

### Element types (`el_types`)

The element types to be used in a model are defined in the `el_types` cell array of strings. The names follow the Abaqus naming convention. Typically, it is easiest to just give a list of all possible element types for the dimensionality of model and then pick the indices of the ones you want to use for each element, e.g.
`el_types = fn_2d_el_types();
el_typ_to_use_for_solid = 'CPE3'; 
mod.el_typ_i(solid_el_indices) = find(strcmp(el_types, el_typ_to_use_for_solid))`
where `solid_el_indices` is a list of the indices of the elements to which you want to assign `CPE3` elements (which are 3-noded plane strain triangular elements).

### Loading steps (`steps`)

The loads that will be applied to a model are defined in the cell array `steps`. Each `step` describes a loading history, in `step{s}.load`, that starts from the original model in its quiescent state - they steps are *not* applied sequentially despite what the name suggests. Each step also defines, in `step{s}.mon`, what will be output ('mon' = monitored) from the solver during that loading.

Typical contents of `step{s}.load`:
- `step{s}.load.time` - `1 x n_time_pts` vector of time steps at which the model will be executed in this loading step
- `step{s}.load.frc_nds` - `n_frc_nds x 1` vector of node indices where loads will be applied (`n_frc_nds` is the number of nodes at which forcing will be applied)
- `step{s}.load.frc_dfs` - `n_frc_nds x 1` vector of the associated Degree of Freedom (DoF) where loads will be applied
- `step{s}.load.frcs` - `1 x n_time_pts` or `n_frc_nds x n_time_pts` matrix or vector of the forcing histories to be applied. If it is a vector, then same forcing history is applied at all nodes/DoFs.
- `steps{1}.load.wts` - `n_frc_nds x 1` optional vector of weightings to be applied to forces at each node/DoF
- 
The requested results for the corresponding loading step are returned in the cell array `res`. Typical outputs are one or both of: 
- History outputs - complete time histories of the displacement (or pressure in fluids) at one or mode nodes, typically plotted as time-domain signals.
- Field output - snapshots of the complete wavefield (its local kinetic energy) at intervals in time, typically displayed as a movie and used as a visualisation tool.

Most of the code in the example scripts is concerned with preparing mod, matls, and steps before fn_FE_entry_point is called and then displaying the outputs.


