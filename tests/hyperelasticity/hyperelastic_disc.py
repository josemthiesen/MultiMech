# Routine to test a hyperelastic disc

from dolfin import *

from mpi4py import MPI

import ufl_legacy as ufl

import numpy as np

import matplotlib.pyplot as plt

from mshr import *

#import periodic_structure as mesher

import source.tool_box.mesh_handling_tools as mesh_tools

import source.constitutive_models.hiperelasticity.isotropic_hyperelasticity as constitutive_models

import source.tool_box.functional_tools as functional_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

import source.tool_box.tensor_tools as tensor_tools

import source.tool_box.variational_tools as variational_tools

########################################################################
########################################################################
##                      User defined parameters                       ##
########################################################################
########################################################################

########################################################################
#                         Material properties                          #
########################################################################

# Sets the Young modulus and the Poisson ratio

E = 100E6

v = 0.4

# Sets the material as a neo-hookean material using the corresponding
# class

constitutive_model = constitutive_models.Neo_Hookean([E, v])

########################################################################
#                                 Mesh                                 #
########################################################################

# Defines the name of the file to save the mesh in. Do not write the fi-
# le termination, e.g. .msh or .xdmf; both options will be saved automa-
# tically

file_name = "tests//test_meshes//disc_mesh"

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

# Sets the solver parameters

linear_solver = "minres"

relative_tolerance = 1e-2

absolute_tolerance = 1e-2

maximum_iterations = 50

preconditioner = "petsc_amg"

krylov_absoluteTolerance = 1e-6

krylov_relativeTolerance = 1e-6

krylov_maximumIterations = 15000

krylov_monitorConvergence = False

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

maximum_load = 2E2

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
boundary_meshCollection, boundary_meshFunction) = mesh_tools.read_xdmfMesh(
file_name)

########################################################################
#                            Function space                            #
########################################################################

# Defines the finite element spaces for the displacement field, u

U = VectorFunctionSpace("CG", mesh.ufl_cell(), polynomial_degree)

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

dsol = TrialFunction(U) 

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

residual_derivative = derivative(residual_form , u_new, dsol)

Res = NonlinearVariationalProblem(residual_form, u_new, bc, J=
residual_derivative)

########################################################################
#                      Solver parameters setting                       #
########################################################################

solver = NonlinearVariationalSolver(Res)

solver.parameters["nonlinear_solver"] = "newton"

solver.parameters["newton_solver"]["linear_solver"] = linear_solver

solver.parameters["newton_solver"]["relative_tolerance"] = (
relative_tolerance)

solver.parameters["newton_solver"]["absolute_tolerance"] = (
absolute_tolerance)

solver.parameters["newton_solver"]["maximum_iterations"] = (
maximum_iterations)

solver.parameters["newton_solver"]["preconditioner"] = (
preconditioner)

solver.parameters['newton_solver']['krylov_solver']['absolute_tole'+
'rance'] = krylov_absoluteTolerance

solver.parameters['newton_solver']['krylov_solver']['relative_tole'+
'rance'] = krylov_relativeTolerance

solver.parameters['newton_solver']['krylov_solver']['maximum_itera'+
'tions'] = krylov_maximumIterations

solver.parameters['newton_solver']['krylov_solver']['monitor_conve'+
'rgence'] = krylov_monitorConvergence

########################################################################
#                         Files initialization                         #
########################################################################

displacement_file = File("./ResultsDir/u.pvd")

########################################################################
#                   Solution and pseudotime stepping                   #
########################################################################

# Initializes the pseudotime counter

time_counter = 0

# Evaluates the pseudotime step

delta_t = (t_final-t)/maximum_loadingSteps

# Iterates through the pseudotime stepping

while t<t_final:

    print("###########################################################"+
    "#############\n#                 Incremental step: "+str(
    time_counter+1)+"; current time: "+str(t)+"               #\n#####"+
    "#################################################################"+
    "##\n")

    # Solves the nonlinear variational problem 

    solver.solve()

    u_new.rename("DNS Displacement", "DNS")

    # Updates the files

    displacement_file << u_new

    # Updates the pseudo time variables, the load, and the counter

    t += delta_t

    load.t = t
    
    time_counter += 1

    if time_counter>=maximum_loadingSteps:

        print("The maximum number of loading steps,",
        maximum_loadingSteps, "has just been reached. Stops the simula"+
        "tion immediatly\n")

        break

print ("\nSimulation completed")