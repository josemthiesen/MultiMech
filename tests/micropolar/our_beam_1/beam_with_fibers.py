# Routine to test a hyperelastic disc

import os

import sys

import numpy as np

from dolfin import *

from mshr import *

import source.constitutive_models.hyperelasticity.micropolar_hyperelasticity as micropolar_constitutiveModels

import source.physics.hyperelastic_micropolar_continuum as variational_framework

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

    # Sets the x, y, and z indices of the RVE to be selected for homoge-
    # nization. These indices begin with 1

    RVE_localizationX = np.ceil(0.5*n_RVEsX)

    RVE_localizationY = np.ceil(0.5*n_RVEsY)

    RVE_localizationZ = np.ceil(0.5*n_RVEsZ)

    # Defines a list of lists, each list is a set of material parameters
    # - micropolar number of the matrix, micropolar number of the fiber,
    # characteristic length of the matrix, characteristic length of the
    # fiber, and load factor

    parameters_sets = [[0.002, 0.002, RVE_width*n_RVEsZ*2.0, (RVE_width*
    n_RVEsZ*2.0), 5.0], [0.02, 0.02, RVE_width*n_RVEsZ*2.0, (RVE_width*
    n_RVEsZ*2.0), 5.0], [0.2, 0.2, RVE_width*n_RVEsZ*2.0, (RVE_width*
    n_RVEsZ*2.0), 5.0]]

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
        ][2], parameters_sets[i][3], True, parameters_sets[i][4], 
        gamma_matrix=0.0, gamma_fiber=0.0, RVE_width=RVE_width, 
        RVE_length=RVE_length, fiber_radius=fiber_radius, n_RVEsX=
        n_RVEsX, n_RVEsY=n_RVEsY, n_RVEsZ=n_RVEsZ, RVE_localizationX=
        RVE_localizationX, RVE_localizationY=RVE_localizationY, 
        RVE_localizationZ=RVE_localizationZ, flag_newMesh=flag_mesh)

# Defines a function to try different parameters

def beam_case_1(E_matrix, E_fiber, nu_matrix, nu_fiber, 
N_micropolarMatrix, N_micropolarFiber, characteristic_lengthMatrix, 
characteristic_lengthFiber, flag_bending, load_factor, gamma_matrix=0.0, 
gamma_fiber=0.0, RVE_width=1.0, RVE_length=1.0, fiber_radius=0.25, 
n_RVEsX=1, n_RVEsY=1, n_RVEsZ=5, RVE_localizationX=1, RVE_localizationY=
1, RVE_localizationZ=3, flag_newMesh=True):
    
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

    displacement_fileName = ["displacement.xdmf", "microrotation.xdmf"]

    homogenized_displacementFileName = ["homogenized_displacement.txt", 
    "homogenized_microrotation.txt"]

    homogenized_gradDisplacementFileName = ["homogenized_displacement_"+
    "gradient.txt", "homogenized_microrotation_grad.txt"]

    stress_fieldFileName = ["cauchy_stress.xdmf", "couple_cauchy_stres"+
    "s.xdmf"]

    post_processes = []

    # Iterates through the fields (displacement and microrotation)

    for i in range(2):

        post_processes.append(dict())

        post_processes[-1]["SaveField"] = {"directory path":
        results_pathGraphics, "file name":displacement_fileName[i]}

        # Put "" in the subdomain to integrate over the entire domain

        post_processes[-1]["HomogenizeField"] = {"directory path":
        results_pathText, "file name":homogenized_displacementFileName[i
        ], "subdomain":""}

        # Put "" in the subdomain to integrate over the entire domain

        post_processes[-1]["HomogenizeFieldsGradient"] = {"directory p"+
        "ath":results_pathText, "file name":
        homogenized_gradDisplacementFileName[i], "subdomain":""}

        # Adds the stress field to the displacement field even though
        # it can be evaluated with any field, since it takes all fields
        # simultaneously

        if i==3:

            post_processes[-1]["SaveStressField"] = {"directory path":
            results_pathText, "file name": stress_fieldFileName[i], 
            "polynomial degree": 1}

    post_processesSubmesh = []

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

    constitutive_model[("Matrix","RVE matrix")] = micropolar_constitutiveModels.Micropolar_Neo_Hookean(
    material_propertiesMatrix)

    constitutive_model[("Fiber","RVE fiber")] = micropolar_constitutiveModels.Micropolar_Neo_Hookean(
    material_propertiesFiber)

    ####################################################################
    #                               Mesh                               #
    ####################################################################

    # Defines the name of the file to save the mesh in. Do not write the 
    # file termination, e.g. .msh or .xdmf; both options will be saved 
    # automatically

    file_directory = os.getcwd()+"//tests//test_meshes"

    mesh_fileName = "micropolar_beam_with_fibers"

    if flag_newMesh:

        beam_gmsh.case_1(RVE_width, RVE_length, fiber_radius, n_RVEsX, 
        n_RVEsY, n_RVEsZ, RVE_localizationX, RVE_localizationY, 
        RVE_localizationZ, mesh_fileName=mesh_fileName, file_directory=
        file_directory)

    # Defines a set of physical groups to create a submesh

    volume_physGroupsSubmesh = []

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

    # Defines a load expression
    
    maximum_load = ((0.5*load_factor*E_matrix*((RVE_length*(RVE_width**3
    ))/12))/((n_RVEsZ*RVE_width)**3))

    maximum_load1 = 0.84E4

    maximum_load = 1E5

    print(maximum_load1/maximum_load)

    load = Expression("(t/t_final)*maximum_load", t=t, t_final=t_final,
    maximum_load=maximum_load, degree=2)

    # Assembles this load into the list of Neumann boundary conditions

    neumann_loads = [load]

    # Assemble the traction vector using this load expression

    traction_boundary = as_vector([0.0, load, 0.0])

    # Defines a dictionary of tractions

    traction_dictionary = dict()

    traction_dictionary["upper"] = traction_boundary

    traction_dictionary["lower"] = traction_boundary

    # Defines a dictionary of moments on the boundary

    moment_boundary = as_vector([0.0, 0.0, 0.0])

    moment_dictionary = dict()

    moment_dictionary["lower"] = moment_boundary

    # Defines the boundary physical groups to apply fixed support boun-
    # dary condition. This variable can be either a list of physical 
    # groups tags or simply a tag. Applies for both displacement and mi-
    # crorotation

    fixed_supportDisplacementPhysicalGroups = "back"

    fixed_supportMicrorotationPhysicalGroups = "back"

    ####################################################################
    ####################################################################
    ##                    Calculation and solution                    ##
    ####################################################################
    ####################################################################

    # Defines a flag to print every step

    verbose = True

    # Solves the variational problem

    variational_framework.hyperelasticity_displacementMicrorotationBased(
    constitutive_model, traction_dictionary, moment_dictionary, 
    maximum_loadingSteps, t_final, post_processes, file_directory+"//"+
    mesh_fileName, solver_parameters, neumann_loads=neumann_loads, 
    polynomial_degreeDisplacement=polynomial_degreeDisplacement, 
    polynomial_degreeMicrorotation=polynomial_degreeMicrorotation,
    t=t, fixed_supportDisplacementPhysicalGroups=
    fixed_supportDisplacementPhysicalGroups, solution_name=[[
    "displacement", "DNS"], ["microrotation", "DNS"]], 
    volume_physGroupsSubmesh=volume_physGroupsSubmesh, 
    fixed_supportMicrorotationPhysicalGroups=
    fixed_supportMicrorotationPhysicalGroups, post_processesSubmesh=
    post_processesSubmesh, verbose=verbose)

case1_varyingMicropolarNumber(flag_newMesh=True)