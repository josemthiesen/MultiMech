# Routine to test a hyperelastic disc

from dolfin import *

import os

from mpi4py import MPI

from mshr import *

#import periodic_structure as mesher

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.file_handling_tools as file_handling_tools

import source.constitutive_models.hiperelasticity.anisotropic_hyperelasticity as anisotropic_constitutiveModels

import source.constitutive_models.hiperelasticity.isotropic_hyperelasticity as isotropic_constitutiveModels

import source.tool_box.functional_tools as functional_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

import source.tool_box.tensor_tools as tensor_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

########################################################################
########################################################################
##                      User defined parameters                       ##
########################################################################
########################################################################

########################################################################
#                          Simulation results                          #
########################################################################

# Defines the path to the results directory 

results_path = os.getcwd()+"//tests//hyperelasticity//results"

displacement_fileName = "displacement.pvd"

########################################################################
#                         Material properties                          #
########################################################################

# Sets a dictionary of properties

material_properties1 = dict()

# Shearing modulus

material_properties1["c"] = 10E6

# k1 is the fiber modulus and k2 is the exponential coefficient

material_properties1["k1"] = 10E4

material_properties1["k2"] = 5.0

# Kappa is the fiber dispersion and it is bounded between 0 and 1/3. A 
# third is an isotropic material

material_properties1["kappa"] = 0.2

# Gamma is the fiber angle in degrees

material_properties1["gamma"] = 45.0

# k is the matrix bulk modulus

material_properties1["k"] = 15E6

# The vectors ahead form a plane where the fiber is locally present

material_properties1["local system of coordinates: a direction"] = (
as_vector([1.0, 0.0, 0.0]))

material_properties1["local system of coordinates: d direction"] = (
as_vector([0.0, 0.0, 1.0]))

material_properties2 = dict()

material_properties2["E"] = 1E7

material_properties2["v"] = 0.4

material_properties3 = dict()

material_properties3["E"] = 1E8

material_properties3["v"] = 0.35

# Sets the material as a HGO material

constitutive_model = dict()

option = 2

if option==1:

    constitutive_model[1] = anisotropic_constitutiveModels.Holzapfel_Gasser_Ogden_Unconstrained(
    material_properties1)

    constitutive_model[2] = isotropic_constitutiveModels.Neo_Hookean(
    material_properties2)

    constitutive_model[3] = isotropic_constitutiveModels.Neo_Hookean(
    material_properties3)

elif option==2:

    constitutive_model[1] = anisotropic_constitutiveModels.Holzapfel_Gasser_Ogden_Unconstrained(
    material_properties1)

    constitutive_model[tuple([2,3])] = isotropic_constitutiveModels.Neo_Hookean(
    material_properties2)

elif option==3:

    constitutive_model = anisotropic_constitutiveModels.Holzapfel_Gasser_Ogden_Unconstrained(
    material_properties1)

########################################################################
#                                 Mesh                                 #
########################################################################

# Defines the name of the file to save the mesh in. Do not write the fi-
# le termination, e.g. .msh or .xdmf; both options will be saved automa-
# tically

file_name = "tests//test_meshes//intervertebral_disc"

########################################################################
#                            Function space                            #
########################################################################

# Defines the shape functions degree

polynomial_degree = 1

########################################################################
#                           Solver parameters                          #
########################################################################

# Sets some parameters

parameters["form_compiler"]["representation"] = "uflacs"
parameters["allow_extrapolation"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 2

# Sets the solver parameters in a dictionary

solver_parameters = dict()

solver_parameters["linear_solver"] = "minres"

solver_parameters["newton_relative_tolerance"] = 1e-4

solver_parameters["newton_absolute_tolerance"] = 1e-4

solver_parameters["newton_maximum_iterations"] = 50

solver_parameters["preconditioner"] = "petsc_amg"

solver_parameters["krylov_absolute_tolerance"] = 1e-6

solver_parameters["krylov_relative_tolerance"] = 1e-6

solver_parameters["krylov_maximum_iterations"] = 15000

solver_parameters["krylov_monitor_convergence"] = False

# Sets the initial time

t = 0.0

# Sets the final pseudotime of the simulation

t_final = 1.0

# Sets the maximum number of steps of loading

maximum_loadingSteps = 30

########################################################################
#                          Boundary conditions                         #
########################################################################

# Defines a load expression

maximum_load = 8E6

load = Expression("(t/t_final)*maximum_load", t=t, t_final=t_final,
maximum_load=maximum_load, degree=0)

# Assembles this load into the list of Neumann boundary conditions

neumann_loads = [load]

# Assemble the traction vector using this load expression

traction_boundary = as_vector([0.0, 0.0, load])

# Defines a dictionary of tractions

traction_dictionary = dict()

traction_dictionary[5] = traction_boundary

# Defines the boundary physical groups to apply fixed support boundary
# condition. This variable can be either a list of physical groups tags
# or simply a tag

fixed_supportPhysicalGroups = 4

########################################################################
########################################################################
##               Calculation: Mr. User, take care ahead!              ##
########################################################################
########################################################################

########################################################################
#                                 Mesh                                 #
########################################################################

# Reads the mesh and constructs some fenics objects using the xdmf file

(mesh, dx, ds, n, domain_meshCollection, domain_meshFunction, 
boundary_meshCollection, boundary_meshFunction) = mesh_tools.read_mshMesh(
file_name)

########################################################################
#                            Function space                            #
########################################################################

# Defines the finite element spaces for the displacement field, u

U = VectorFunctionSpace(mesh, "Lagrange", polynomial_degree)

########################################################################
#                          Boundary conditions                         #
########################################################################

# Defines the boundary conditions for fixed facets

bc = BCs_tools.fixed_supportDirichletBC(U, boundary_meshFunction, 
fixed_supportPhysicalGroups)

########################################################################
#                           Variational forms                          #
########################################################################

# Defines the trial and test functions

delta_u = TrialFunction(U) 

v = TestFunction(U)

# Creates the function for the updated solution, i.e. the vector of pa-
# rameters

u_new = Function(U)

# Constructs the variational form for the inner work

internal_VarForm = variational_tools.hyperelastic_internalWorkFirstPiola(
u_new, v, constitutive_model, dx)

# Constructs the variational forms for the traction work

traction_VarForm = variational_tools.traction_work(traction_dictionary,
v, ds)

# Assembles the residual, takes the Gateaux derivative and assembles the
# nonlinear problem object

residual_form = internal_VarForm-traction_VarForm

residual_derivative = derivative(residual_form , u_new, delta_u)

Res = NonlinearVariationalProblem(residual_form, u_new, bc, J=
residual_derivative)

########################################################################
#                      Solver parameters setting                       #
########################################################################

solver = NonlinearVariationalSolver(Res)

########################################################################
#                         Files initialization                         #
########################################################################

# Creates the path to the displacement file

displacement_file = file_handling_tools.verify_path(results_path, 
displacement_fileName)

displacement_file = File(displacement_file)

########################################################################
#                   Solution and pseudotime stepping                   #
########################################################################

# Evaluates the pseudotime step

delta_t = (t_final-t)/maximum_loadingSteps

# Iterates through the pseudotime stepping algortihm 

newton_raphson_tools.newton_raphsonSingleField(t, t_final, delta_t, 
maximum_loadingSteps, solver, u_new, displacement_file, neumann_loads=
neumann_loads, solver_parameters=solver_parameters)