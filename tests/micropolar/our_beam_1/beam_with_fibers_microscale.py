# Routine to test a hyperelastic disc

import os

import sys

import numpy as np

from dolfin import *

from mshr import *

import source.constitutive_models.hyperelasticity.micropolar_hyperelasticity as micropolar_constitutiveModels

import source.multiscale.multiscale_micropolar as variational_framework

from source.tool_box.file_handling_tools import float_toString

sys.path.insert(1, '/home/matheus-janczkowski/Github')

import CuboidGmsh.tests.micropolar_meshes.beam_micropolar_case_1 as beam_gmsh

# Defines a function to try multiple parameters

def case1_varyingMicropolarNumber(flag_newMesh=False):

    # Sets the Young modulus and the Poisson ration from psi to MPa

    E_matrix = 100E6

    nu_matrix = 0.4

    E_fiber = 1000E6

    nu_fiber = 0.4

    # Defines the RVE overall parameters

    RVE_width = 1.0

    RVE_length = 1.0

    # Defines the fiber radius

    fiber_radius = 0.25

    # Defines the number of RVEs at each direction

    n_RVEsX = 1

    n_RVEsY = 1

    n_RVEsZ = 5

    # Defines a list of lists, each list is a set of material parameters
    # - micropolar number of the matrix, micropolar number of the fiber,
    # characteristic length of the matrix, characteristic length of the
    # fiber, and load factor

    parameters_sets = [[0.002, 0.002, RVE_width*n_RVEsZ*2.0, (RVE_width*
    n_RVEsZ*2.0), 5.0]]#, [0.02, 0.02, RVE_width*n_RVEsZ*2.0, (RVE_width*
    #n_RVEsZ*2.0), 5.0]]#, [0.2, 0.2, RVE_width*n_RVEsZ*2.0, (RVE_width*
    #n_RVEsZ*2.0), 5.0]]

    # Iterates through the simulations

    for i in range(len(parameters_sets)):

        # Makes a new mesh just for the first test and if a new mesh is
        # asked for

        flag_mesh = False

        if flag_newMesh and i==0:

            flag_mesh = True

        # Calls the simulation for bending

        beam_case_1(E_matrix, E_fiber, nu_matrix, nu_fiber, 
        parameters_sets[i][0], parameters_sets[i][1], parameters_sets[i
        ][2], parameters_sets[i][3], True, gamma_matrix=0.0, gamma_fiber
        =0.0, RVE_width=RVE_width, RVE_length=RVE_length, fiber_radius=
        fiber_radius, flag_newMesh=flag_mesh)

# Defines a function to try different parameters

def beam_case_1(E_matrix, E_fiber, nu_matrix, nu_fiber, 
N_micropolarMatrix, N_micropolarFiber, characteristic_lengthMatrix, 
characteristic_lengthFiber, flag_bending, gamma_matrix=0.0, gamma_fiber=
0.0, RVE_width=1.0, RVE_length=1.0, fiber_radius=0.25, n_RVEsX=1, 
n_RVEsY=1, n_RVEsZ=1, RVE_localizationX=1, RVE_localizationY=1, 
RVE_localizationZ=3, flag_newMesh=True):
    
    # Sets a sufix to denote the material parameters

    sufix = ""

    if flag_bending:

        sufix += ("bending_matrix_lb_"+float_toString(
        characteristic_lengthMatrix)+"_fiber_lb_"+float_toString(
        characteristic_lengthFiber)+"_")

    else:

        sufix += ("torsion_matrix_lt_"+float_toString(
        characteristic_lengthMatrix)+"_fiber_lt_"+float_toString(
        characteristic_lengthFiber)+"_")

    sufix += ("matrix_gamma_"+float_toString(gamma_matrix)+"_fiber_gam"+
    "ma_"+float_toString(gamma_fiber)+"_")

    sufix += ("matrix_N_"+float_toString(N_micropolarMatrix)+"_fiber_N"+
    "_"+float_toString(N_micropolarFiber))

    ####################################################################
    ####################################################################
    ##                    User defined parameters                     ##
    ####################################################################
    ####################################################################

    ####################################################################
    #                        Simulation results                        #
    ####################################################################

    # Defines the path to the results directory 

    results_pathGraphics = (os.getcwd()+"//tests//micropolar//our_beam"+
    "_1//results//graphics//"+sufix)

    results_pathText = (os.getcwd()+"//tests//micropolar//our_beam_1//"+
    "results//text//"+sufix)

    displacement_fileName = ["displacement_microscale.xdmf", "microrot"+
    "ation_microscale.xdmf"]

    homogenized_displacementFileName = ["homogenized_displacement_micr"+
    "oscale.txt", "homogenized_microrotation_microscale.txt"]

    homogenized_gradDisplacementFileName = ["homogenized_displacement_"+
    "gradient_microscale.txt", "homogenized_microrotation_grad_microsc"+
    "ale.txt"]

    stress_fieldFileName = ["cauchy_stress_microscale.xdmf", "couple_c"+
    "auchy_stress_microscale.xdmf"]

    post_processes = []

    # Iterates through the fields (displacement and microrotation)

    for i in range(2):

        post_processes.append([i, dict()])

        post_processes[-1][-1]["SaveField"] = {"directory path":
        results_pathGraphics, "file name":displacement_fileName[i]}

        # Put "" in the subdomain to integrate over the entire domain

        post_processes[-1][-1]["HomogenizeField"] = {"directory path":
        results_pathText, "file name":homogenized_displacementFileName[i
        ], "subdomain":""}

        # Put "" in the subdomain to integrate over the entire domain

        post_processes[-1][-1]["HomogenizeFieldsGradient"] = {"directo"+
        "ry path":results_pathText, "file name":
        homogenized_gradDisplacementFileName[i], "subdomain":""}

        # Adds the stress field to the displacement field even though
        # it can be evaluated with any field, since it takes all fields
        # simultaneously

        if i==0:

            post_processes[-1][-1]["SaveStressField"] = {"directory pa"+
            "th": results_pathGraphics, "file name": 
            stress_fieldFileName[0], "polynomial degree": 1}

            post_processes[-1][-1]["SaveCoupleStressField"] = {"direct"+
            "ory path": results_pathGraphics, "file name": 
            stress_fieldFileName[1], "polynomial degree": 1}

    ####################################################################
    #                       Material properties                        #
    ####################################################################

    # Converts to Lamé parameters

    mu_matrix = E_matrix/(2*(1+nu_matrix))

    lmbda_matrix = (nu_matrix*E_matrix)/((1+nu_matrix)*(1-(2*nu_matrix)))

    mu_fiber = E_fiber/(2*(1+nu_fiber))

    lmbda_fiber = (nu_fiber*E_fiber)/((1+nu_fiber)*(1-(2*nu_fiber)))

    # Sets a dictionary of properties

    alpha_matrix = 0.0

    beta_matrix = 0.0

    alpha_fiber = 0.0

    beta_fiber = 0.0

    if flag_bending:

        beta_matrix = 4*mu_matrix*(characteristic_lengthMatrix**2)

        beta_fiber = 4*mu_fiber*(characteristic_lengthFiber**2)

    else:

        beta_matrix = ((2*mu_matrix*(characteristic_lengthMatrix**2))-
        gamma_matrix)

        beta_fiber = ((2*mu_fiber*(characteristic_lengthFiber**2))-
        gamma_fiber)

    # Saves the properties into a dictionary for the matrix

    material_propertiesMatrix = dict()

    material_propertiesMatrix["mu"] = mu_matrix

    material_propertiesMatrix["lambda"] = lmbda_matrix

    material_propertiesMatrix["N"] = N_micropolarMatrix

    material_propertiesMatrix["alpha"] = alpha_matrix

    material_propertiesMatrix["beta"] = beta_matrix

    material_propertiesMatrix["gamma"] = gamma_matrix

    # And for the fiber

    material_propertiesFiber = dict()

    material_propertiesFiber["mu"] = mu_fiber

    material_propertiesFiber["lambda"] = lmbda_fiber

    material_propertiesFiber["N"] = N_micropolarFiber

    material_propertiesFiber["alpha"] = alpha_fiber

    material_propertiesFiber["beta"] = beta_fiber

    material_propertiesFiber["gamma"] = gamma_fiber

    # Sets the material as a HGO material

    constitutive_model = dict()

    constitutive_model["Matrix"] = micropolar_constitutiveModels.Micropolar_Neo_Hookean(
    material_propertiesMatrix)

    constitutive_model["Fiber"] = micropolar_constitutiveModels.Micropolar_Neo_Hookean(
    material_propertiesFiber)

    ####################################################################
    #                               Mesh                               #
    ####################################################################

    # Defines the name of the file to save the mesh in. Do not write the 
    # file termination, e.g. .msh or .xdmf; both options will be saved 
    # automatically

    file_directory = os.getcwd()+"//tests//test_meshes"

    mesh_fileName = "micropolar_beam_with_fibers_microscale"

    if flag_newMesh:

        beam_gmsh.case_1(RVE_width, RVE_length, fiber_radius, n_RVEsX, 
        n_RVEsY, n_RVEsZ, RVE_localizationX, RVE_localizationY, 
        RVE_localizationZ, mesh_fileName=mesh_fileName, file_directory=
        file_directory)

    ####################################################################
    #                          Function space                          #
    ####################################################################

    # Defines the shape functions degree

    polynomial_degreeDisplacement = 2

    polynomial_degreeMicrorotation = 1

    ####################################################################
    #                         Solver parameters                        #
    ####################################################################

    # Sets the solver parameters in a dictionary

    solver_parameters = dict()

    solver_parameters["linear_solver"] = "mumps"#"mumps"

    solver_parameters["newton_relative_tolerance"] = 1e-8#1e-8

    solver_parameters["newton_absolute_tolerance"] = 1e-5#1e-8

    solver_parameters["newton_maximum_iterations"] = 30

    """

    solver_parameters["preconditioner"] = "petsc_amg"

    solver_parameters["krylov_absolute_tolerance"] = 1e-5

    solver_parameters["krylov_relative_tolerance"] = 1e-5

    solver_parameters["krylov_maximum_iterations"] = 15000

    solver_parameters["krylov_monitor_convergence"] = False"""

    # Sets the initial time

    t = 0.0

    # Sets the final pseudotime of the simulation

    t_final = 1.0

    # Sets the maximum number of steps of loading

    maximum_loadingSteps = 25

    ####################################################################
    #                        Boundary conditions                       #
    ####################################################################

    # Defines the paths to the macro quantities files

    macro_displacementName = (results_pathText+"//homogenized_displace"+
    "ment")

    macro_gradDisplacementName= (results_pathText+"//homogenized_displ"+
    "acement_gradient") 

    macro_microrotationName = (results_pathText+"//homogenized_microro"+
    "tation") 

    macro_gradMicrorotationName = (results_pathText+"//homogenized_mic"+
    "rorotation_grad")

    ####################################################################
    ####################################################################
    ##                    Calculation and solution                    ##
    ####################################################################
    ####################################################################

    # Defines a flag to print every step

    verbose = True

    # Solves the variational problem

    variational_framework.micropolar_microscale(macro_displacementName, 
    macro_gradDisplacementName, macro_microrotationName, 
    macro_gradMicrorotationName, constitutive_model,
    maximum_loadingSteps, t_final, post_processes, file_directory+"//"+
    mesh_fileName, solver_parameters, polynomial_degreeDisplacement=
    polynomial_degreeDisplacement, polynomial_degreeMicrorotation=
    polynomial_degreeMicrorotation, t=t, solution_name=[["displacement", 
    "DNS"], ["microrotation", "DNS"]], verbose=verbose)

case1_varyingMicropolarNumber(flag_newMesh=True)