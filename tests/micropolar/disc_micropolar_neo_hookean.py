# Routine to test a hyperelastic disc

from dolfin import *

import os

#from mpi4py import MPI

from mshr import *

#import periodic_structure as mesher

import source.constitutive_models.hyperelasticity.micropolar_hyperelasticity as micropolar_constitutiveModels

import source.physics.hyperelastic_micropolar_continuum as variational_framework

########################################################################
########################################################################
##                      User defined parameters                       ##
########################################################################
########################################################################

########################################################################
#                          Simulation results                          #
########################################################################

# Defines the path to the results directory 

results_path = os.getcwd()+"//tests//micropolar//results"

displacement_fileName = ["displacement.xdmf", "microrotation.xdmf"]

homogenized_displacementFileName = ["homogenized_displacement.txt", "h"+
"omogenized_microrotation.txt"]

homogenized_gradDisplacementFileName = ["homogenized_displacement_grad"+
"ient.txt", "homogenized_microrotation_grad.txt"]

post_processes = []

# Iterates through the fields (displacement and microrotation)

for i in range(2):

    post_processes.append(dict())

    post_processes[-1]["save field"] = {"directory path":results_path, 
    "file name":displacement_fileName[i]}

    # Put "" in the subdomain to integrate over the entire domain

    post_processes[-1]["homogenize field"] = {"directory path":
    results_path, "file name":homogenized_displacementFileName[i], 
    "subdomain":""}

    # Put "" in the subdomain to integrate over the entire domain

    post_processes[-1]["homogenize field's gradient"] = {"directory path":
    results_path, "file name":homogenized_gradDisplacementFileName[i], 
    "subdomain":""}

post_processesSubmesh = []

########################################################################
#                         Material properties                          #
########################################################################

# Sets a dictionary of properties

material_properties = dict()

mu = 26.12

material_properties["mu"] = mu

material_properties["lambda"] = 63.84-(2*mu/3)

material_properties["kappa"] = 0.0

material_properties["alpha"] = 0.0

material_properties["beta"] = 0.0

material_properties["gamma"] = 1e-12

# Sets the material as a HGO material

constitutive_model = dict()

constitutive_model[tuple([1,2,3])] = micropolar_constitutiveModels.Micropolar_Neo_Hookean(
material_properties)

########################################################################
#                                 Mesh                                 #
########################################################################

# Defines the name of the file to save the mesh in. Do not write the fi-
# le termination, e.g. .msh or .xdmf; both options will be saved automa-
# tically

mesh_fileName = "tests//test_meshes//intervertebral_disc"

# Defines a set of physical groups to create a submesh

volume_physGroupsSubmesh = []

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

solver_parameters["linear_solver"] = "mumps"#"minres"

solver_parameters["newton_relative_tolerance"] = 1e-10#1e-3

solver_parameters["newton_absolute_tolerance"] = 1e-10#1e-3

solver_parameters["newton_maximum_iterations"] = 10#50

"""

solver_parameters["preconditioner"] = "petsc_amg"

solver_parameters["krylov_absolute_tolerance"] = 1e-6

solver_parameters["krylov_relative_tolerance"] = 1e-6

solver_parameters["krylov_maximum_iterations"] = 15000

solver_parameters["krylov_monitor_convergence"] = False"""

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

maximum_load = 4E1

load = Expression("(t/t_final)*maximum_load", t=t, t_final=t_final,
maximum_load=maximum_load, degree=0)

# Assembles this load into the list of Neumann boundary conditions

neumann_loads = [load]

# Assemble the traction vector using this load expression

traction_boundary = as_vector([0.0, 0.0, load])

# Defines a dictionary of tractions

traction_dictionary = dict()

traction_dictionary[5] = traction_boundary

# Defines a dictionary of moments on the boundary

moment_boundary = as_vector([0.0, 0.0, 0.0])

moment_dictionary = dict()

moment_dictionary[5] = moment_boundary

# Defines the boundary physical groups to apply fixed support boundary
# condition. This variable can be either a list of physical groups tags
# or simply a tag. Applies for both displacement and microrotation

fixed_supportDisplacementPhysicalGroups = 4

fixed_supportMicrorotationPhysicalGroups = 4

########################################################################
########################################################################
##                      Calculation and solution                      ##
########################################################################
########################################################################

# Solves the variational problem

variational_framework.hyperelasticity_displacementMicrorotationBased(
constitutive_model, traction_dictionary, moment_dictionary, 
maximum_loadingSteps, t_final, post_processes, mesh_fileName, 
solver_parameters, neumann_loads=neumann_loads, polynomial_degree=
polynomial_degree, t=t, fixed_supportDisplacementPhysicalGroups=
fixed_supportDisplacementPhysicalGroups, solution_name=[["displacement",
"DNS"], ["microrotation", "DNS"]], volume_physGroupsSubmesh=
volume_physGroupsSubmesh, fixed_supportMicrorotationPhysicalGroups=
fixed_supportMicrorotationPhysicalGroups, post_processesSubmesh=
post_processesSubmesh)