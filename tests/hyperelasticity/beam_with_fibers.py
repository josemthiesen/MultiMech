# Routine to test a hyperelastic disc

from dolfin import *

import os

from mpi4py import MPI

from mshr import *

#import periodic_structure as mesher

import source.constitutive_models.hyperelasticity.isotropic_hyperelasticity as constitutive_models

import source.physics.hyperelastic_cauchy_continuum as variational_framework

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

displacement_fileName = "displacement.xdmf"

post_processes = dict()

post_processes["save field"] = {"directory path":results_path, 
"file name":displacement_fileName}

########################################################################
#                         Material properties                          #
########################################################################

# Sets the Young modulus and the Poisson ratio

E = 100E6

nu = 0.4

# Sets a dictionary of properties

material_properties = dict()

material_properties["E"] = E

material_properties["v"] = nu

# Sets the material as a neo-hookean material using the corresponding
# class

constitutive_model = constitutive_models.Neo_Hookean(material_properties)

########################################################################
#                                 Mesh                                 #
########################################################################

# Defines the name of the file to save the mesh in. Do not write the fi-
# le termination, e.g. .msh or .xdmf; both options will be saved automa-
# tically

mesh_fileName = "tests//test_meshes//micropolar_beam_with_fibers"

########################################################################
#                            Function space                            #
########################################################################

# Defines the shape functions degree

polynomial_degree = 1

########################################################################
#                           Solver parameters                          #
########################################################################

# Sets the solver parameters in a dictionary

solver_parameters = dict()

solver_parameters["linear_solver"] = "mumps"

solver_parameters["newton_relative_tolerance"] = 1e-4

solver_parameters["newton_absolute_tolerance"] = 1e-4

solver_parameters["newton_maximum_iterations"] = 10

"""solver_parameters["preconditioner"] = "hypre_amg"

solver_parameters["krylov_absolute_tolerance"] = 1e-9

solver_parameters["krylov_relative_tolerance"] = 1e-9

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

maximum_load = -2E6#1E-1

load = Expression("(t/t_final)*maximum_load", t=0.0, t_final=t_final,
maximum_load=maximum_load, degree=0)

# Assembles this load into the list of Neumann boundary conditions

neumann_loads = [load]

# Assemble the traction vector using this load expression

traction_boundary = as_vector([0.0, load, 0.0])

# Defines a dictionary of tractions

traction_dictionary = dict()

traction_dictionary["bottom"] = traction_boundary

# Defines a load expression for prescribed displacement

maximum_displacement = 1E-1

displacement_load = Expression("(t/t_final)*maximum_displacement", t=t,
t_final=t_final, maximum_displacement=maximum_displacement, degree=0)

# Assembles this prescribed displacement into the list of Dirichlet 
# boundary conditions

dirichlet_loads = [displacement_load]

# Sets the dictionary of prescribed displacement. Each key is a physical
# group or a tuple of physical groups, and the value is a list in either
# one of the following formats:
# value = [load_expression]        -> applies the load to all DOFs
# value = [1, load_expression]     -> applies the load to one DOF
# value = [[0,1], load_expression] -> applies the load to a list of DOFs

prescribed_displacement = dict()

#prescribed_displacement["bottom"] = [1, displacement_load]

# Defines the boundary physical groups to apply fixed support boundary
# condition. This variable can be either a list of physical groups tags
# or simply a tag

fixed_supportPhysicalGroups = "back"

########################################################################
########################################################################
##                      Calculation and solution                      ##
########################################################################
########################################################################

# Solves the variational problem

variational_framework.hyperelasticity_displacementBased(
constitutive_model, traction_dictionary, maximum_loadingSteps, t_final, 
post_processes, mesh_fileName, solver_parameters, neumann_loads=
neumann_loads, dirichlet_loads=dirichlet_loads, prescribed_displacement=
prescribed_displacement, polynomial_degree=polynomial_degree, t=t, 
fixed_supportPhysicalGroups=fixed_supportPhysicalGroups, verbose=True)