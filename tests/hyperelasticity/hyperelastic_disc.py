# Routine to test a hyperelastic disc

from dolfin import *

import os

from mpi4py import MPI

import ufl_legacy as ufl

import numpy as np

import matplotlib.pyplot as plt

from mshr import *

#import periodic_structure as mesher

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.file_handling_tools as file_handling_tools

import source.constitutive_models.hiperelasticity.isotropic_hyperelasticity as constitutive_models

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

# Sets the Young modulus and the Poisson ratio

E = 100E6

v = 0.4

# Sets a dictionary of properties

material_properties = dict()

material_properties["E"] = E

material_properties["v"] = v

# Sets the material as a neo-hookean material using the corresponding
# class

constitutive_model = constitutive_models.Neo_Hookean(material_properties)

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

maximum_loadingSteps = 11

########################################################################
#                          Boundary conditions                         #
########################################################################

# Defines a load expression

maximum_load = 2E7

load = Expression("(t/t_final)*maximum_load", t=t, t_final=t_final,
maximum_load=maximum_load, degree=0)

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

# Creates the deformation gradient object

I = Identity(3)

F = grad(u_new)+I

# Initializes objects for the stresses at the reference configuration

first_piola = constitutive_model.first_piolaStress(F)

# Constructs the variational forms for the inner work

internal_VarForm = inner(first_piola, grad(v))*dx 

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
maximum_loadingSteps, solver, u_new, displacement_file, neumann_loads=[
load], solver_parameters=solver_parameters)