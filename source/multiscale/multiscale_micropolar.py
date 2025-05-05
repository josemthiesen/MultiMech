# Routine to store the multiscale methods for micropolar physics

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

import source.tool_box.programming_tools as programming_tools

# Defines a function to model a hyperelastic problem with a displacement
# and a microrotation fields only in the microscale. It uses the macro
# quantities read from txt files

@programming_tools.optional_argumentsInitializer({'neumann_loads': 
lambda: [], 'dirichlet_loads': lambda: [], 'solution_name': lambda: [], 
'simple_supportDisplacementPhysicalGroups': lambda: dict(), ('simple_s'+
'upportMicrorotationPhysicalGroups'): lambda: dict(), ('volume_physGro'+
'upsSubmesh'): lambda: [], 'post_processesSubmesh': lambda: []})

def micropolar_microscale(macro_displacementName, 
macro_gradDisplacementName, macro_microrotationName, 
macro_gradMicrorotationName, constitutive_model, traction_dictionary, 
moment_dictionary, maximum_loadingSteps, t_final, post_processes, 
mesh_fileName, solver_parameters, polynomial_degreeDisplacement=2, 
polynomial_degreeMicrorotation=2, t=0.0, 
fixed_supportDisplacementPhysicalGroups=0, neumann_loads=None, 
dirichlet_loads=None, fixed_supportMicrorotationPhysicalGroups=0, 
solution_name=None, simple_supportDisplacementPhysicalGroups=None, 
simple_supportMicrorotationPhysicalGroups=None, volume_physGroupsSubmesh
=None, post_processesSubmesh=None, verbose=False):

    ####################################################################
    #                               Mesh                               #
    ####################################################################

    # Reads the mesh and constructs some fenics objects using the xdmf 
    # file

    mesh_dataClass = mesh_tools.read_mshMesh(mesh_fileName, verbose=
    verbose)

    ####################################################################
    #                          Function space                          #
    ####################################################################

    # Constructs elements for the displacement and for the microrotation
    # fields

    displacement_element = VectorElement("CG", 
    mesh_dataClass.mesh.ufl_cell(), polynomial_degreeDisplacement)

    microrotation_element = VectorElement("CG", 
    mesh_dataClass.mesh.ufl_cell(), polynomial_degreeMicrorotation)

    # Defines the finite element spaces for the Lagrange multipliers, 
    # which enforce the micro-scale kinematical constraints. The degree 
    # of the polynomial is 0 because the lagrange multiplier is constant
    # throughout the element

    # Displacement constraint (impedes rigid body translation)

    lagrangeMult_U = VectorElement("Real", mesh_dataClass.mesh.ufl_cell(
    ), 0) 

    # Displacement gradient constraint (impedes rigid body rotation)

    lagrangeMult_GradU = TensorElement("Real", 
    mesh_dataClass.mesh.ufl_cell(), 0) 

    # Microrotation constraint

    lagrangeMult_Phi = VectorElement("Real", 
    mesh_dataClass.mesh.ufl_cell(), 0) 

    # Microrotation gradient constraint

    lagrangeMult_GradPhi = TensorElement("Real", 
    mesh_dataClass.mesh.ufl_cell(), 0)

    # Defines the finite element spaces for post-processing

    Space_Piola_1st = TensorElement("DG", mesh_dataClass.mesh.ufl_cell(
    ), 0)

    # Defines the mixed element for the monolithic solution

    mixed_element = MixedElement([displacement_element, 
    microrotation_element, lagrangeMult_U, lagrangeMult_GradU, 
    lagrangeMult_Phi, lagrangeMult_GradPhi])

    # Defines the finite element space for the monolithic solution

    monolithic_functionSpace = FunctionSpace(mesh_dataClass.mesh, 
    mixed_element)

    ####################################################################
    #                         Macro quantities                         #
    ####################################################################

    # Reads the macro displacement, microrotation, and their gradients

    u_RVE_homogenized = file_tools.txt_toList(macro_displacementName)

    grad_u_RVE_homogenized = file_tools.txt_toList(
    macro_gradDisplacementName)

    phi_RVE_homogenized = file_tools.txt_toList(macro_microrotationName)

    grad_phi_RVE_homogenized = file_tools.txt_toList(
    macro_gradMicrorotationName)

    ####################################################################
    #                         Variational forms                        #
    ####################################################################

    # Defines the trial and test functions

    delta_solution = TrialFunction(monolithic_functionSpace) 

    variation_solution = TestFunction(monolithic_functionSpace)

    # Creates the function for the updated solution, i.e. the vector of 
    # parameters

    solution_new = Function(monolithic_functionSpace)

    # Splits the solution and the test function. Splits the fields but 
    # keeps each one with the global vector of parameters

    u_new, phi_new = split(solution_new)

    variation_u, variation_phi = split(variation_solution)

    # Constructs the variational form for the inner work

    internal_VarForm = variational_tools.hyperelastic_micropolarInternalWorkFirstPiola(
    u_new, phi_new, variation_u, variation_phi, constitutive_model, 
    mesh_dataClass)

    # Constructs the variational forms for the traction work

    traction_VarForm = variational_tools.traction_work(
    traction_dictionary, variation_u, mesh_dataClass)

    # Constructs the variational forms for the moment work on the boun-
    # dary. Note that the function traction_work was reused, because the
    # variational construction is the same for traction and for moment

    moment_VarForm = variational_tools.traction_work(
    moment_dictionary, variation_phi, mesh_dataClass)

    # Assembles the residual, takes the Gateaux derivative and assembles
    # the nonlinear problem object

    residual_form = internal_VarForm-traction_VarForm-moment_VarForm

    residual_derivative = derivative(residual_form, solution_new, 
    delta_solution)

    # Makes the boundary condition an empty list, because the macrosca-
    # le boundary conditions are applied directly onto the variational
    # form, using Lagrange multipliers

    Res = NonlinearVariationalProblem(residual_form, solution_new, [], 
    J=residual_derivative)

    ####################################################################
    #                    Solver parameters setting                     #
    ####################################################################

    solver = NonlinearVariationalSolver(Res)

    ####################################################################
    #                 Solution and pseudotime stepping                 #
    ####################################################################

    # Iterates through the pseudotime stepping algortihm 

    newton_raphson_tools.newton_raphsonMultipleFields(t, t_final, 
    maximum_loadingSteps, solver, solution_new, mixed_element, 
    mesh_dataClass, constitutive_model, post_processesList=
    post_processes, post_processesSubmeshList=post_processesSubmesh, 
    dirichlet_loads=dirichlet_loads, neumann_loads=neumann_loads, 
    solver_parameters=solver_parameters, volume_physGroupsSubmesh=
    volume_physGroupsSubmesh, solution_name=solution_name)