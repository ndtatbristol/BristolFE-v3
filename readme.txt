BRISTOLFE v3 - Paul Wilcox
==========================

This repository contains a number of Matlab functions and example scripts for performing basic Finite Element (FE) simulations, in particular explicit ones for elastodynamic wave propagation. Simulations can be performed using the built-in solver or a Pogo solver if installed (NB Pogo does not currently support a pressure formulation for fluids, hence its use is limited to solid-only simulations). An interface to Abaqus may be added at some point. The built-in solver is currently limited to 2D models are and the following types of element: CPE3 3-noded triangular elements for elastic solids, AC2D3 3-noded triangular elements for fluids, ASI2D2 2-noded elements for fluid-solid interface elements.

To use, clone (or download and unzip) the repository and add the 'BristolFE-v3/code' folder to the Matlab path.

The entry-point for FE analysis is the function 'fn_FE_entry_point' in the 'BristolFE-v3/code' folder. The 'BristolFE-v2/code' folder also contains numerous helper functions for creating meshes, adding defects, and displaying results.

The core FE code is in various functions in the 'BristolFE-v2/code/internal' folder, which are not expected to be called directly by the user. 

Special wrapper functions for sub-domain modelling are in the 'BristolFE-v3/subdoms' folder. See the two example files for examples of how to use.

Make sure that the code and code/internal folders are on the Matlab path, e.g. by having:
addpath(genpath('RELEVANT_PATH/BristolFE-v2/code'));
at the top of any scripts that use the functions, where RELEVANT_PATH is set according to wherever you put the folder.

Some example scripts are provided in the 'BristolFE-v3/examples'. Most likely you will start with a copy of one of these and modify it according to your requirements.

SUMMARY OF CHANGES SINCE v2
===========================

The v2 to v3 changes are designed to make the code more consistent and future-proof for expansion to 3D. The main changes are:

- The entry point for all solvers is now 'fn_FE_entry_point' (not 'fn_BristolVE_v2'), with the solver being selected as an option.
- 'matls' is now a cell array rather than a structure array. Each element in the model (except for interface elements) remains associated with a material via the n_elsx1 vector, 'mod.el_mat_i', which contains integers referencing the corresponding cell in 'matls'.
- 'el_types' is a new cell array of strings, one for each type of element used in the model. The element index vector, , that associates each element with an element type is now a vector of integer references
- 'mod.el_typ_i' is no longer a cell array of strings describing element types, but instead n_elsx1 vector of integer indices that reference the corresponding cell in 'el_types'.
- 'mod.el_abs_i' is still an n_elsx1 vector of relative absorption of elements in the range 0 to 1 but is now mandatory in the mod structure, not optional
- Names and locations of functions in 'BristolFE-v3/code' are being rationalised, and are in the process of having consistent help comments added. All functions not intended for direct use in the 'BristolFE-v3/code/internal folder'. Functions exclusively for use in 2D models are prefixed 'fn_2d_'; functions exclusively for use in 3D models are prefixed by 'fn_2d_'. Functions without these prefixes can be used in 2D or 3D models with the same arguments.

To facilitate transition to v3, modified versions of deprecated v2 functions that take the v2 input and output arguments are included in 'BristolFE-v3/code/deprecated'; if called, these issue a warning, internally modify the input parameters if necessary, call the appropriate v3 function, modify that function's output if necessary, and return the same output as the original v2 function. This means that v2 scripts should still run without modification and for this reason the v2 examples in their original form are included in 'BristolFE-v3/v2 examples'. A new set of equivalent examples are now provided in 'BristolFE-v2/examples' that call the appropriate v3 functions directly and should be used as the basis for future scripts.

OVERVIEW
========

The entry point function is res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options). 

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


