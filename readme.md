# BristolFE v3.1

### Paul Wilcox

This repository contains a number of Matlab functions and example scripts for performing explicit time-marching Finite Element (FE) simulations for elastodynamic wave propagation. Simulations can be performed using the built-in solver that exploits Matlab's native sparse-matrix and GPU processing capabilities or an alternative commercial solver such as [*pogo*](www.pogo.software). Support for Abaqus may be added in the future. 

## Installation

To use, clone (or download and unzip) the repository and either permanently add the `BristolFE-v3/code` folder to the Matlab path or include the line

```
addpath(genpath('RELEVANT_PATH/BristolFE-v3/code'));
```

at the start of any scripts that use the *BristolFE* functions.

## Getting started

The user interacts with *BristolFE* via the suite of documented Matlab functions in the folder `BristolFE-v3/code`. These provide tools for the creation, execution, and visualisation of models. They are intended to be called from user-written Matlab script files and a number of example scripts are provided in `BristolFE-v3/examples`. Most likely you will start with a copy of one of these and modify it according to your requirements.

## Examples

In `BristolFE-v3/examples` you will find examples how to set up different features in Bristol FE models. The file `very_basic_example_in_single_script.m` is a good starting point as it is a self-contained script that demonstrates how to create, execute, and visualise the results from a model with minimum complexity. However, it is generally better practice to define a model parametrically in a suitable function to keep code tidy. Example model functions, `mod_*.m`, can be found in `BristolFE-v3/examples/models`. The script `general_examples.m` can be used to run any of these and visualise the results - just uncomment the line for the relevant model file. The general format of the example model functions is `[mod, matls, el_types, steps, fe_options, params] = mod_NAME(params)`, where `params` is a structured variable containing the parameters of a given model. The returned values, `mod`, `matls`, `el_types`, `steps`, and `fe_options` are ready to be used in a call to `fn_FE_entry_point` function. Each model file has a complete set of default parameters (`default_params`) specified at the top of the file, which will be used unless they are overwritten by any of the fields in the `params` argument passed to the model function. For example, setting `params.els_per_wavelength = 6` (elements per wavelength) will overwrite the `default_params.els_per_wavelength` value specified in the model function. The other returned value, `params`, from a `mod_*.m` function contains a complete list of the parameters used for the model (i.e. the default parameters plus any overwritten by user). The example model files include:
- `mod_2d_basic.m` - a parameterised version of the simple model in `very_basic_example_in_single_script.m`
- `mod_2d_advanced.m` - a 2D model which can include fluid and solid domains, absorbing boundaries, volumetric scatterers, and a zero-volume crack, all of which can be turned on or off by setting the appropriate params.
- `mod_2d_oblique.m`- a classic NDE configuration of a transducer on a perspex wedge designed to insonify surface-breaking cracks on the far side of a steel plate.
- `mod_3d_basic.m` - a simple 3D solid model
- `mod_3d_advanced.m` - a more complex 3D model, with roughly similar features to `mod_2d_advanced.m`
- `mod_2d_subdomain.m` - a 2D model that demonstrates the concept of sub-domains (run using the `subdomain_examples.m` script)
- `mod_2d_subdomain_oblique.m` - a 2D sub-domain model version of `mod_2d_oblique.m` (run using the `subdomain_examples.m` script)

## Summary of changes since v3.0

### Element library extended

Element files are now built in a self-consistent way that covers 2D, 3D, and interface elements. The following elements are available:
- `CPE3` - 3-noded 2D solid, linear, triangular
- `CPE4` - 4-noded 2D solid, linear, quadrilateral
- `AC2D3` - 3-noded 2D fluid, linear, triangular
- `AC2D4` - 4-noded 2D fluid, linear, quadrilateral
- `ASID2` - 2-noded 2D fluid-solid interface
- `C3D4` - 4-noded 3D solid, linear, tetrahedral
- `C3D6` - 6-noded 3D solid, linear, triangular prism
- `C3D8` - 8-noded 3D solid, linear, hexahedral
- `AC3D4` - 4-noded 3D fluid, linear, tetrahedral
- `AC3D8` - 8-noded 3D fluid, linear, hexahedral
- `ASI3D3` - 3-noded 3D fluid-solid interface, triangular
- `ASI3D4` - 4-noded 3D fluid-solid interface, quadrilateral

### Addition of a material library function

To save time defining common materials, a utility function, `matl = fn_material_library(name)`, has been added that returns a material structure to be used in models based on the specified name, which must be one of the materials in the library. The list will be expanded in future versions.

### Improved global matrix builder (v7)

This has been re-written to improve efficiency and avoid out-of-memory errors in large models. Out-of-memory errors will still occur for very large models, but the previous version caused them to occur unnecessarily during the matrix building phase.

### Solver options in main time-marching routine (v7)

A predictor-corrector solver mode is now available which may be more stable than the standard explicit scheme in models with fluid-solid coupling. The available solver modes are:
- `explicit` / `exp` / `vel at last half time step` - is fastest but may be unstable in models with fluid-solid coupling unless the time-step is very short
- `implicit` / `imp` / `vel at current time step` - is guaranteed to be stable but requires initial inversion of non-diagonal damping matrix which may be prohibitively slow, especially in 3D models
- `predictor corrector` / `pc` - is as fast as `explicit` and should be a bit more stable in coupled fluid-solid problems

## Summary of changes since v2

The v2 to v3 changes are designed to make the code more consistent and future-proof for expansion to 3D. The main changes are:
- The built-in solver has been made more efficient and the previous issue of instability when using fluid-solid interface elements has been resolved by reformulating the solver to use velocity at the current time step rather than at the previous half-step
- Names of functions have been rationalised: functions exclusively for use in 2D models are prefixed `fn_2d_`; functions exclusively for use in 3D models are prefixed by `fn_3d_`; functions without these prefixes can be used in 2D or 3D models with the same arguments
- The entry point for all solvers is now `fn_FE_entry_point` (not `fn_BristolVE_v2`), with the solver being selected as an option (the default is `fe_options.solver = 'BristolFE'`).
- A cell array, `matls`, is now used to describe materials rather than a structure array. Each element in the model (except for interface elements) remains associated with a material via the `n_els`×1 vector, `mod.el_mat_i`, which contains integers referencing the corresponding cell in `matls`
- `el_types` is a new cell array of strings, one for each type of element used in the model
- `mod.el_typ_i` is no longer a cell array of strings describing element types, but instead `n_els`×1 vector of integer indices that reference the corresponding cell in `el_types`.
- `mod.el_abs_i` is still an`n_els`×1 vector of relative absorption of elements in the range 0 (no absorption) to 1 (maximum absorption) but is now mandatory in the `mod` structure, not optional

## Overview

The entry point function for solving a model is the function `fn_FE_entry_point`, which typically appears in a script like this:

```
%Code to prepare the model before this point

%Execute the model
res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%Code to process the results from the model after this point
``` 

When `fn_FE_entry_point` is called, a complete mesh must have been specified (in the structure `mod`), the materials used must have been defined (in the cell array `matls`), the element types used must have been defined (in the cell array `el_types`) and one or more loading steps and the required outputs must have been defined (in the cell array `steps`). 

### Execution options (`fe_options`)

The parameter `fe_options` is a structure that can contain the following fields (if a field is not part of the structure, the indicated default is used):
- `fe_options.solver = [default = 'BristolFE']` - which solver to use. Current options are `'BristolFE'` or `'pogo'`
- `fe_options.solver_precision [default = 'double']` - precision used for calculations.  Options are `'single'` or `'double'`
- `fe_options.field_output_every_n_frames [default = inf]` - specifies how often field output should be recorded, if at all. See below.
- `fe_options.dof_to_use [default = []]` - which Degrees-of-Freedom (DoF) to include in calculations. Leave empty for all. 
- `fe_options.sort_nds [default = 0]` - set to 1 to sort order of nodes based on physical coordinates prior to solving. Seems to be needed for some models executed in *pogo*.
- `fe_options.nd_sort_cols [default = []]` - if `fe_options.sort_nds` = 1, this sets the column order of coordinates if they are being sorted.
- `fe_options.damping_power_law [default = 3]`  - the power law index that governs the rate of damping increase in a absorbing layer as the absorption index goes from 0 to 1
- `fe_options.max_damping [default = []]` - the maximum damping to be applied in a absorbing layer (i.e. at an absorption index of 1). Leave empty for default with is the reciprocal of the time-step used
- `fe_options.max_stiffness_reduction [default = 0.01]` - the maxium reduction in stiffness to be applied in a absorbing layer, i.e. at an absorption index of 1

Options specific to *BristolFE* solver:

- `fe_options.use_gpu_if_available [default = 1]` - as it says, use the GPU if there is one available
- `fe_options.solver_mode [default = 'implicit']` - determines how nodal velocity is calculated in the time-marching solver. The default `'implicit'` is always stable, but less efficient than `'explicit'` which often leads to instability if a fluid-solid boundary is present in the model
- `fe_options.field_output_type [default = 'KE']` - the nodal quantity output whenever a frame of field output is recorded. Currently the only option is `'KE'` (kinetic energy)

Options specific to the *pogo* solver:

- `fe_options.pogo_path [default = 'C:\Program Files\pogo\windows']` - location of the *pogo* solver executable
- `fe_options.pogo_matlab_path [default = 'C:\Program Files\pogo\matlab']` - location of the *pogo* Matlab functions
- `fe_options.pogo_verbosity [default = -1]` - level of *pogo* output to Matlab console (-1 = none)
- `fe_options.pogo_compression [default = 0]` - whether to allow *pogo* to do compression by approximating values for similar elements
- `fe_options.pogo_number_of_diff_absorbing_matls [default = inf]` - *pogo* uses different materials with the approriate damping levels to create an absorbing layer (whereas BristolFE adds damping to the elements in an absorbing layer based on their absorbing index at runtime). This option determines how many different absorbing materials *pogo* should use for absorbing boundary layers. If the default of `'inf'` is used then a different material will be created for every unique value of absorbing index found, which should gave identical result to *BristolFE*

### Model description (`mod`)

This describes the model geometry and must contain the following fields:
- `mod.nds` - an `n_nds`×`n_dims` matrix of coordinates of all nodes in the model (`n_nds` is number of nodes; `n_dims` is number of dimensions, i.e. 2 or 3)
- `mod.els` - an `n_els`×`max_nds_per_el` matrix of the nodes associated with each element in the model (`n_els` is number of elements; `max_nds_per_el` is the maximum number of nodes used by an element in the model, which is typically 3 for a 2D model of triangular elements but can be as many as 8 for a 3D model with hexahedral elements)
- `mod.el_mat_i` - an `n_els`×1 vector of indices that describe the material associated with each element
- `mod.el_typ_i` - an `n_els`×1 vector of indices that describe the element type associated with each element
- `mod.el_abs_i` - an `n_els`×1 vector of numbers in the range 0 to 1 that describe the relative damping of each element (damping can be specified directly as a material property, but the relative damping controls how the elements damping and stiffness behaviour is modified when it is part of an absorbing layer)

*Note that node and element numbers are implicitly defined by the associated row number of the relevant matrix, where the first row represents node or element 1 (not 0).*
 
 ### Material descriptions (`matls`)
 The materials to be used in a model are defined in a `matls` cell array, with the cell index being the identifier references by `mod.el_mat_i`. The requried fields are:
 - `matls{i}.name` - a string giving the name
 - `matls{i}.rho` - the density
 - `matls{i}.D` - either:
   - solids - a 6×6 stiffness matrix (in a 2D model the necessary reduction to a 3×3 stiffness matrix takes place when plane stress or plane straing elements are formed)
   - fluids - a 1×1 matrix containing the bulk modulus
 - `matls{i}.col` - an 1×3 vector of RBG values used to determine the colour used for plotting. 

Functions such as the following are provided that convert engineering data into the necessary material property values for some common cases:
- `matls{i} = fn_matl_isotropic_solid_defined_by_velocities(name, longitudinal_velocity, shear_velocity, density)`
- `matls{i} = fn_matl_fluid_defined_by_velocity(name, velocity, density)`

### Element types (`el_types`)

The element types to be used in a model are defined in the `el_types` cell array of strings. The names follow the Abaqus naming convention. Typically, it is easiest to just give a list of all possible element types for the dimensionality of model and then pick the indices of the ones you want to use for each element, e.g.
```
el_types = fn_2d_el_types();
el_typ_to_use_for_solid = 'CPE3'; 
mod.el_typ_i(solid_el_indices) = find(strcmp(el_types, el_typ_to_use_for_solid))`
```
where `solid_el_indices` is a list of the indices of the elements to which you want to assign `CPE3` elements (which are 3-noded plane strain triangular elements).

### Loading steps (`steps`)

The loads that will be applied to a model are defined in the cell array `steps`. Each `step` describes a loading history, in `step{s}.load`, that starts from the original model in its quiescent state - they steps are *not* applied sequentially despite what the name suggests. Each step also defines, in `step{s}.mon` ('mon' = monitored), what history values will be recorded from the solver during that loading. Each loading and history monitoring point is defined by a node/Degree-of-Freedom (DoFs) pair. At this point, it is necessary to introduce the DoFs used in *BristolFE*, which are:
1. force or displacement in the x-direction in solids
2. force or displacement in the y-direction in solids
3. force or displacement in the z-direction in solids
4. volumetric expansion rate or time-integral of pressure (proportional to velocity potential) in fluids

In the above list, the first term for each DoF is the type of load that can be applied to the DoF and the second term is the type of output that can be recorded for that DoF. (Internally, the latter are the field quantities that the FE solver solves for at every node/DoF pair throughout the model at every time point.)  

Typical contents of `step{s}.load`:
- `step{s}.load.time` - a 1×`n_time_pts` vector of time steps at which the model will be executed in this loading step
- `step{s}.load.frc_nds` - an `n_frc_nds`×1 vector of node indices where loads will be applied (`n_frc_nds` is the number of nodes at which forcing will be applied)
- `step{s}.load.frc_dfs` - an `n_frc_nds`×1 vector of the associated DoFs where loads will be applied
- `step{s}.load.frcs` - a 1×`n_time_pts` or `n_frc_nds`×`n_time_pts` matrix or vector of the forcing histories to be applied. If it is a vector, then same forcing history is applied at all nodes/DoFs.
- `steps{s}.load.frc_wts` - an `n_frc_nds`×1 optional vector of weightings to be applied to forces at each node/DoF. This provides an efficient way of applying a single load that is not aligned to a single DoF direction at each node while still only requiring a vector for `steps{s}.load.frcs` (note that `wts` is now deprecated and `frc_wts` should be used instead although `wts` will still work for the moment)

#### Requesting history outputs

The contents of `step{s}.mon` define what history outputs will be recorded in the output by specifying the node/DoF pairs in a similar manner to the way loads are specified:
- `step{s}.mon.dsp_nds` - an `n_mon_nds`×1 vector of node indices where history values will be recorded (`n_mon_nds` is the number of nodes at which history outputs are requested) (note that `nds` is now deprecated and `dsp_nds` should be used instead although `nds` will still work for the moment)
- `step{s}.mon.dsp_dfs` - an `n_mon_nds`×1 vector of the associated DoFs for which history values will be recorded (note that `dfs` is now deprecated and `dsp_dfs` should be used instead although `dfs` will still work for the moment)
- `step{s}.mon.dsp_wts` - an optional `n_output`×`n_mon_nds` matrix of weightings to get linear combination(s) of the monitored nodal displacements rather than those values directly, effectively complementing the concept of frc_wts for loading

#### Requesting field output

To record field outputs (snapshots of the field values at all nodes in the model), set `fe_options.field_output_every_n_frames` to a small integer. For example, if set to 10, a snapshot will be recorded every 10 time-steps. Obviously, a smaller number will result in more data, but around something in the 5-20 range is usually fine. If you don't want field ouptuts, use `fe_options.field_output_every_n_frames = inf`. This is the default if not specified. 

### Model results (`res`)

The output from a model is stored in the cell array `res`, in which each cell corresponds to the associated loading step. Typical contents of a `res{s}`:
- `res{s}.dsps` - an `n_frc_nds`×`n_time_pts` matrix of the requested history outputs
- `res{s}.valid_mon_dsps` - an `n_frc_nds`×1 logical vector with 1 wherever a requested history output is valid (which should be all of them)
- `res{s}.fld` - an `n_total_dofs`×`n_fld_time_pts` matrix of field outputs if requested (`n_total_dofs` is the total number of DoFs in the entire model, `n_fld_time_pts` are the number of field outputs (determined by value of `fe_options.field_output_every_n_frames`)
- `res{s}.fld_time`  - a 1×`n_fld_time_pts` vector of times associated with each field output

Typically, field outputs are used for visualisation and are used in conjuntion with the *BristolFE* functions `fn_show_geometry` and `fn_run_animation` in code like:

```
figure;
display_options.draw_elements = 0;
h_patch = fn_show_geometry(mod, matls, display_options);
anim_options.repeat_n_times = 1;
fn_run_animation(h_patch, res{1}.fld, anim_options);
```
## Function summary
### fn_2d_add_absorbing_layer

Adds an absorbing boundary by increasing element absorbing indices proportional to their distance from the specified boundary divided by the specified absorbing boundary layer thickness (i.e. so it reaches one when the distance is equal to the absorbing boundary layer thickness. The boundary defines the start of the absorbing layer; within the boundary the absorbing index is set to zero.
```
mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness [, els_to_apply_to])
```
### fn_2d_add_crack

Adds a crack into a 2D by identifying nearest element edges/faces and 'splitting' model along them, by duplicating nodes. Default is a zero width crack unless optional Crack Opening Displacement (COD) is specified in which case the nodes are displaced away from plane of crack
```
mod = fn_2d_add_crack(mod, el_types, crack_vtcs, [crack_fcs, cod])
```
### fn_2d_add_inclusion_or_void

Adds scatterer to existing model by turning all elements inside scat_pts to either matl(scat_matl) or void if = scat_matl
```
mod = fn_2d_add_inclusion_or_void(mod, el_types, scat_pts, scat_matl, scat_el_typ)
```
### fn_2d_create_smooth_random_blob

Creates list of 2D points that describe perimeter of smooth blob shape.
```
pts = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts)
```
### fn_2d_el_types

Returns cell array of available 2D element types. Use in 2D models, usually exactly like this: el_types = fn_2d_el_types()
```
el_types = fn_2d_el_types()
```
### fn_2d_find_elements_in_region

Returns logical n_els x 1 vectors indicating whether elements in model are inside or outside the specified region.
```
[in, out] = fn_2d_find_elements_in_region(mod, region)
```
### fn_2d_random_walk

Returns
```
pts = fn_2d_random_walk(npts, step_mean, step_std, angle_start, angle_mean, angle_std)
```
### fn_2d_rough_line

Returns rough line relative to mean straight line described by rms roughness and correlation length
```
pts = fn_2d_rough_line(length, angle, rms_roughness, corr_len)
```
### fn_2d_signed_dist_to_bdry

Returns signed (positive exterior) shortest distance of point(s) to boundary surface described by vertices of triangular facets
```
[d, nearest_pts, norm_vecs, type_of_nearest_entity, nearest_entity, bdry_edges] = fn_2d_signed_dist_to_bdry(pts, bdry_nds, bdry_edges)
```
### fn_2d_structured_mesh

Utility function for generating a structured mesh of triangular or quadrilateral elements, that fills the region specified by bdry_nds.
```
mod = fn_2d_structured_mesh(bdry_pts, el_size, el_type [, force_equilateral_els])
```
### fn_2d_structured_mesh_rectangular_els

Utility function for generating a structured mesh of square elements
```
mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size)
```
### fn_2d_structured_mesh_triangular_els

Utility function for generating a isometric structured mesh of triangular elements, that fills the region specified by bdry_nds.
```
mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size [, force_equilateral_els])
```
### fn_3d_add_absorbing_layer

Adds an absorbing boundary by increasing element absorbing indices proportional to their distance from the specified boundary divided by the specified absorbing boundary layer thickness (i.e. so it reaches one when the distance is equal to the absorbing boundary layer thickness. The boundary defines the start of the absorbing layer; within the boundary the absorbing index is set to zero.
```
mod = fn_3d_add_absorbing_layer(mod, abs_bdry_nds, abs_bdry_fcs, abs_bdry_thickness [, els_to_apply_to])
```
### fn_3d_add_crack

Adds a crack into a 3D model by identifying nearest element edges/faces and 'splitting' model along them, by duplicating nodes. Default is a zero width crack unless optional Crack Opening Displacement (COD) is specified in which case the nodes are displaced away from plane of crack
```
mod = fn_3d_add_crack(mod, el_types, crack_vtcs, crack_fcs [, cod])
```
### fn_3d_add_inclusion_or_void

Adds scatterer to existing model by turning all elements inside scat_pts to either matl(scat_matl) or void if = scat_matl
```
mod = fn_3d_add_inclusion_or_void(mod, el_types, scat_pts, scat_matl, scat_el_typ)
```
### fn_3d_cylindrical_surface

Create a 3D cylindrical surface described by vertices and triangular faces
```
[vtcs, fcs] = fn_3d_cylindrical_surface(pt1, pt2, rad [, n_ang_divisions])
```
### fn_3d_disk_surface

Create a 3D disk surface described by vertices and triangular faces
```
[vtcs, fcs] = fn_3d_disk_surface(pt1, pt2, rad [, n_ang_divisions])
```
### fn_3d_el_types

Returns cell array of available 3D element types. Use in 3D models, usually exactly like this: el_types = fn_3d_el_types()
```
el_types = fn_3d_el_types()
```
### fn_3d_elliptical_surface

Create a 3D elliptical surface described by vertices and triangular faces
```
[vtcs, fcs] = fn_3d_elliptical_surface(pt1, pt2, rad [, n_ang_divisions])
```
### fn_3d_find_elements_in_region

Returns logical n_els x 1 vectors indicating whether elements in model are inside or outside the specified region, which is defined in terms of vertices and triangular faces.
```
[in, out] = fn_3d_find_elements_in_region(mod, vtcs, fcs)
```
### fn_3d_hexahedral_surface

Create a 3D hexahedral surface described by vertices and triangular faces based list of 8 corner vertices, ordered cyclically around 'bottom' face and the cyclically around 'top' face in same sense (i.e. like Abaqus node numbering)
```
[vtcs, fcs] = fn_3d_hexahedral_surface(pts)
```
### fn_3d_rectalinear_surface

Create a 3D rectalinear surface described by vertices and triangular faces based on two opposing corner positions
```
[vtcs, fcs] = fn_3d_rectalinear_surface(pt1, pt2)
```
### fn_3d_signed_dist_to_bdry

Returns signed (positive exterior) shortest distance of point(s) to boundary surface described by vertices of triangular facets
```
d = fn_3d_signed_dist_to_bdry(pts, bdry_nds, bdry_fcs)
```
### fn_3d_spherical_surface

Create a 3D sphere described by vertices and faces
```
[vtcs, fcs] = fn_3d_spherical_surface(cent, rad [, n_sub_divisions])
```
### fn_3d_structured_mesh_hexahedral_els

Utility function for generating a 3d structured mesh of cuboidal elements, that fills the cuboidal region specified by crnr_pts.
```
mod = fn_3d_structured_mesh_hexahedral_els(crnr_pts, el_size)
```
### fn_FE_entry_point

Common entry point for different FE solvers.
```
res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options) [res, mats] = fn_FE_entry_point(mod, matls, el_types, steps, fe_options) fe_options = fn_FE_entry_point([], [], [], fe_options)
```
### fn_add_fluid_solid_interface_els

Adds the necessary interface elements between all solid and fluid elements in a model. Without these there is no coupling between the solid and fluid domains. Only needs to be called once for a given model after all features are added (i.e. it should be the last step)
```
mod = fn_add_fluid_solid_interface_els(mod, el_types)
```
### fn_display_result_v2

Displays mesh from 2D or 3D model, returning handle to patches for later animations
```
fn_display_result(nodes, elements, display_options) to display mesh OR
```
### fn_estimate_matl_vels

Estimates minimum and maximum speeds of sound possible for given material. For isotropic elastic solids these correspond to shear and longitudinal wave speeds. For fluids they will be same and correspond to bulk pressure wave speed. For anisotropic solids, who knows?
```
[max_velocity, min_velocity] = fn_estimate_matl_vels(matls)
```
### fn_estimate_max_min_vels

Estimates minimum and maximum speeds of sound possible given cell array of materials
```
[max_velocity, min_velocity] = fn_estimate_max_min_vels(matls)
```
### fn_estimate_max_min_wavelengths

Estimates minimum and maximum wavelengths possible given cell array of materials
```
[max_wavelength, min_wavelength] = fn_estimate_max_min_wavelengths(matls, centre_frequency)
```
### fn_find_node_nearest_to_point

Finds nearest node to specified point
```
node = fn_find_node_nearest_to_point(nodes, p, tol)
```
### fn_find_nodes_nearest_to_line

Returns list of nodes that lie along line (with specified tolerance) defined by its endpoints
```
[node_list, s, r] = fn_find_nodes_nearest_to_line(nodes, p1, p2, tol)
```
### fn_fluid_el_types

Returns cell array of available fluid element types.
```
el_types = fn_fluid_el_types()
```
### fn_gaussian_pulse

Returns a Gaussian pulse with specifed centre frequency and number of cycles given a specified time axis.
```
s = fn_gaussian_pulse(t, centre_freq, no_cycles[, db_down_at_start, db_down])
```
### fn_get_min_max_element_sizes

Returns minimum and maximum side lengths of elements in a model
```
[min_el_size, max_el_size] = fn_get_min_max_element_sizes(mod, el_types)
```
### fn_get_suitable_el_size

Estimated element size to achieve the desired number of elements per wavelength for the slowest wave possible in any of the materials given. For isotropic solids, the slowest wave will be the shear mode, so for solid models, this function returns an element size based on the number of elements per SHEAR wavelength at the centre frequency.
```
el_size = fn_get_suitable_el_size(matls, nominal_cent_freq, els_per_wavelength)
```
### fn_get_suitable_time_step

Returns what should be a stable time step by calculating the fastest possible wavespeed in the materials, workout out how fast such a wave traverses the specified element size and dividing that by a safety factor (default = sqrt(2), specify a alternative value as 3rd optional argument if desired
```
time_step = fn_get_suitable_time_step(matls, el_size [, safety_factor])
```
### fn_hann_pulse

Returns a Hann-windowed pulse with specifed centre frequency and number of cycles given a specified time axis.
```
s = fn_hann_pulse(t, centre_freq, no_cycles)
```
### fn_hilbert

Performs Hilbert transform
```
h = fn_hilbert(s)
```
### fn_material_library

Returns material structure for material with given name if found in libary
```
matl = fn_material_library(name)
```
### fn_matl_fluid_defined_by_bulk_modulus

Returns material structure for fluid based on specifed bulk_modulus and density.
```
matl = fn_matl_fluid_defined_by_bulk_modulus(name, bulk_modulus, density [, col, options])
```
### fn_matl_fluid_defined_by_velocity

Returns material structure for fluid based on specifed velocity and density.
```
matl = fn_matl_fluid_defined_by_velocity(name, velocity, density [, col, options])
```
### fn_matl_isotropic_solid_defined_by_stiffness

Returns material structure based on specifed youngs_modulus, poissons_ratio, and density.
```
matl = fn_matl_isotropic_solid_defined_by_stiffness(name, youngs_modulus, poissons_ratio, density [, col, options])
```
### fn_matl_isotropic_solid_defined_by_velocities

Returns material structure based on specifed velocities and density.
```
matl = fn_matl_isotropic_solid_defined_by_velocities(name, longitudinal_velocity, shear_velocity, density, [, col, options])
```
### fn_optimum_absorbing_bdry_properties

Returns what appear to be good (not necessarily optimum!) values for the parameters of absorbing boundary layers
```
[max_damping, damping_power_law, max_stiffness_reduction] = fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, matls, centre_freq)
```
### fn_parse_fn_file

Parses function file and extracts contents of each part of help
```
[usage, summary, inputs, outputs] = fn_parse_fn_file(fname)
```
### fn_run_animation

Animates field output on a previously displayed model geometry. The wave intensity is presented as the local kinetic energy on a dB scale (default range is 40dB, normalised to peak amplitude over all locations and all times). The range and normalisation value can be over-ridden by setting the appropriate options (see below).
```
fn_run_animation(h_patch, fld, anim_options)
```
### fn_show_geometry

Plots the model geometry and returns handle to patches used for each element that can be used as an input argument for subsequent field animations if desired.
```
h_patch = fn_show_geometry(mod, matls, el_types, options)
```
### fn_solid_el_types

Returns cell array of available solid element types.
```

```
