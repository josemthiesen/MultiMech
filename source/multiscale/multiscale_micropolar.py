# Routine to store the multiscale methods for micropolar physics

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.functional_tools as functional_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

import source.tool_box.programming_tools as programming_tools

# Defines a function to model a hyperelastic problem with a displacement
# and a microrotation fields only in the microscale. It uses the macro
# quantities read from txt files

@programming_tools.optional_argumentsInitializer({'neumann_loads': 
lambda: [], 'dirichlet_loads': lambda: [], 'solution_name': lambda: []})

def micropolar_microscale(macro_displacementName, 
macro_gradDisplacementName, macro_microrotationName, 
macro_gradMicrorotationName, constitutive_model, maximum_loadingSteps, 
t_final, post_processes, mesh_fileName, solver_parameters, 
polynomial_degreeDisplacement=2, polynomial_degreeMicrorotation=2, t=
0.0, solution_name=None, verbose=False):

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

    # Creates a dictionary with the names of the macro variables (name
    # them as python variables, i.e. no spaces nor non-ASCII characters)
    # as keys and the txt file name

    macro_quantitiesDict = dict()

    macro_quantitiesDict["macro_displacement"] = macro_displacementName

    macro_quantitiesDict["macro_gradDisplacement"] = macro_gradDisplacementName

    macro_quantitiesDict["macro_rotation"] = macro_microrotationName

    macro_quantitiesDict["macro_gradMicrorotation"] =  macro_gradMicrorotationName

    # Initializes the class of macro quantities. Put this class into a
    # list because you could have multiple classes

    macro_quantitiesClasses = [functional_tools.MacroQuantitiesInTime(
    macro_quantitiesDict)]

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

    (u_new, phi_new, lagrange_displacement, lagrange_gradDisp, 
    lagrange_microrotation, lagrange_gradMicrorotation) = split(
    solution_new)

    (variation_u, variation_phi, v_lagrangeDisplacement, v_lagrangeGradDisp, 
    v_lagrangeMicrorotation, v_lagrangeGradMicrorotation) = split(
    variation_solution)

    # Constructs the variational form for the inner work

    internal_VarForm = variational_tools.hyperelastic_micropolarInternalWorkFirstPiola(
    u_new, phi_new, variation_u, variation_phi, constitutive_model, 
    mesh_dataClass)

    # Evaluates the volume of the domain
    
    total_volume = assemble(1.0*mesh_dataClass.dx)

    # Constructs the variational form of the Lagrange multipliers

    Int_dlambu = ((1/total_volume)*((dot(v_lagrangeDisplacement,
    macro_quantitiesClasses[0].macro_displacement)*dx)-(dot(
    v_lagrangeDisplacement,u_new)*dx)))

    Int_dlambphi = ((1/total_volume)*((dot(v_lagrangeMicrorotation,
    macro_quantitiesClasses[0].macro_rotation)*dx)-(dot(
    v_lagrangeMicrorotation, phi_new)*dx)))

    Int_dlambgrau = ((1/total_volume)*((inner(v_lagrangeGradDisp, 
    macro_quantitiesClasses[0].macro_gradDisplacement)*dx)-(inner(
    v_lagrangeGradDisp, grad(u_new))*dx)))

    Int_dlambgradphi = ((1/total_volume)*((inner(
    v_lagrangeGradMicrorotation, macro_quantitiesClasses[0
    ].macro_gradMicrorotation)*dx)-(inner(v_lagrangeGradMicrorotation,
    grad(phi_new))*dx))) 

    # Assembles the residual, takes the Gateaux derivative and assembles
    # the nonlinear problem object

    residual_form = (internal_VarForm-Int_dlambu-Int_dlambphi-
    Int_dlambgrau-Int_dlambgradphi)

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
    post_processes, solver_parameters=solver_parameters, solution_name=
    solution_name, macro_quantitiesClasses=macro_quantitiesClasses)