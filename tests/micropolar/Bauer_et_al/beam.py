# Routine to test a hyperelastic disc

import os

from dolfin import *

from mshr import *

import numpy as np

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

results_pathGraphics = (os.getcwd()+"//tests//micropolar//Bauer_et_al/"+
"/results//graphics")

results_pathText = (os.getcwd()+"//tests//micropolar//Bauer_et_al//res"+
"ults//text")

displacement_fileName = ["displacement.xdmf", "microrotation.xdmf"]

homogenized_displacementFileName = ["homogenized_displacement.txt", "h"+
"omogenized_microrotation.txt"]

homogenized_gradDisplacementFileName = ["homogenized_displacement_grad"+
"ient.txt", "homogenized_microrotation_grad.txt"]

post_processes = []

# Iterates through the fields (displacement and microrotation)

for i in range(2):

    post_processes.append(dict())

    post_processes[-1]["save field"] = {"directory path":results_pathGraphics, 
    "file name":displacement_fileName[i]}

    # Put "" in the subdomain to integrate over the entire domain

    post_processes[-1]["homogenize field"] = {"directory path":
    results_pathText, "file name":homogenized_displacementFileName[i], 
    "subdomain":""}

    # Put "" in the subdomain to integrate over the entire domain

    post_processes[-1]["homogenize field's gradient"] = {"directory path":
    results_pathText, "file name":homogenized_gradDisplacementFileName[i], 
    "subdomain":""}

post_processesSubmesh = []

########################################################################
#                         Material properties                          #
########################################################################

# Sets a dictionary of properties

material_properties = dict()

mu = 26.12

K_constitutive = 63.84

alpha = 0.0

beta = 0.0

beam_sectionWidth = 5E-3

characteristic_length = beam_sectionWidth*(1.5E-1)

kappa = -(0.15)*0.0*mu

# Evaluates gamma

gamma = (2*mu*(characteristic_length**2))-beta

# Saves the properties into a dictionary

material_properties["mu"] = mu

material_properties["lambda"] = K_constitutive-(2*mu/3)

material_properties["kappa"] = kappa

material_properties["alpha"] = alpha

material_properties["beta"] = beta

material_properties["gamma"] = gamma

# Sets the material as a HGO material

constitutive_model = dict()

constitutive_model["Generic volume"] = micropolar_constitutiveModels.Micropolar_Neo_Hookean(
material_properties)

#constitutive_model = micropolar_constitutiveModels.Micropolar_Neo_Hookean(
#material_properties)

########################################################################
#                                 Mesh                                 #
########################################################################

# Defines the name of the file to save the mesh in. Do not write the fi-
# le termination, e.g. .msh or .xdmf; both options will be saved automa-
# tically

mesh_fileName = "tests//test_meshes//micropolar_prism"

# Defines a set of physical groups to create a submesh

volume_physGroupsSubmesh = []

########################################################################
#                            Function space                            #
########################################################################

# Defines the shape functions degree

polynomial_degreeDisplacement = 2

polynomial_degreeMicrorotation = 1

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

solver_parameters["newton_relative_tolerance"] = 1e-8#1e-3

solver_parameters["newton_absolute_tolerance"] = 1e-8#1e-3

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

maximum_load = 4E-2

load = Expression("(t/t_final)*maximum_load", t=t, t_final=t_final,
maximum_load=maximum_load, degree=0)

# Assembles this load into the list of Neumann boundary conditions

neumann_loads = [load]

# Assemble the traction vector using this load expression

traction_boundary = as_vector([load, 0.0, 0.0])

# Defines a dictionary of tractions

traction_dictionary = dict()

traction_dictionary["lower"] = traction_boundary

# Defines a dictionary of moments on the boundary

moment_boundary = as_vector([0.0, 0.0, 0.0])

moment_dictionary = dict()

moment_dictionary["lower"] = moment_boundary

# Defines the boundary physical groups to apply fixed support boundary
# condition. This variable can be either a list of physical groups tags
# or simply a tag. Applies for both displacement and microrotation

fixed_supportDisplacementPhysicalGroups = "back"

fixed_supportMicrorotationPhysicalGroups = "back"

########################################################################
########################################################################
##                      Calculation and solution                      ##
########################################################################
########################################################################

# Solves the variational problem

variational_framework.hyperelasticity_displacementMicrorotationBased(
constitutive_model, traction_dictionary, moment_dictionary, 
maximum_loadingSteps, t_final, post_processes, mesh_fileName, 
solver_parameters, neumann_loads=neumann_loads, 
polynomial_degreeDisplacement=polynomial_degreeDisplacement, 
polynomial_degreeMicrorotation=polynomial_degreeMicrorotation,
t=t, fixed_supportDisplacementPhysicalGroups=
fixed_supportDisplacementPhysicalGroups, solution_name=[["displacement",
"DNS"], ["microrotation", "DNS"]], volume_physGroupsSubmesh=
volume_physGroupsSubmesh, fixed_supportMicrorotationPhysicalGroups=
fixed_supportMicrorotationPhysicalGroups, post_processesSubmesh=
post_processesSubmesh)