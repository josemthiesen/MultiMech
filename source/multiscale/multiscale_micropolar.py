# Routine to store the multiscale methods for micropolar physics

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.functional_tools as functional_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

import source.tool_box.programming_tools as programming_tools

import source.tool_box.multiscale_boundary_conditions_tools as multiscale_BCsTools

# Defines a function to model a hyperelastic problem with a displacement
# and a microrotation fields only in the microscale. It uses the macro
# quantities read from txt files

@programming_tools.optional_argumentsInitializer({'neumann_loads': 
lambda: [], 'dirichlet_loads': lambda: [], 'solution_name': lambda: []})

def micropolar_microscale(displacement_multiscaleBC, 
microrotation_multiscaleBC, macro_displacementFileName, 
macro_gradDisplacementFileName, macro_microrotationFileName, 
macro_gradMicrorotationFileName, constitutive_model, maximum_loadingSteps, 
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

    # Assembles a dictionary of finite elements for the two primal 
    # fields: displacement and microrotation

    elements_dictionary = {"displacement": VectorElement("CG", 
    mesh_dataClass.mesh.ufl_cell(), polynomial_degreeDisplacement), "m"+
    "icrorotation": VectorElement("CG", mesh_dataClass.mesh.ufl_cell(), 
    polynomial_degreeMicrorotation)}

    ####################################################################
    #                        Boundary conditions                       #
    ####################################################################

    # Builds a dictionary of multiscale boundary conditions. Each key is
    # the name of the field to be constrained, whereas the value is a
    # dictionary with two mandatory keys: "boundary condition" that sta-
    # tes the type of boundary condition to be applied and has the name
    # of one of the classes defined in multiscale_classes.py; "macro in-
    # formation", that is a dictionary of the macro file paths

    multiscale_BCsDict = dict()

    multiscale_BCsDict["displacement"] = {"boundary condition": 
    displacement_multiscaleBC, "macro information": {"macro field file":
    macro_displacementFileName, "macro field gradient file":
    macro_gradDisplacementFileName}}

    multiscale_BCsDict["microrotation"] = {"boundary condition": 
    microrotation_multiscaleBC, "macro information": {"macro field fil"+
    "e": macro_microrotationFileName, "macro field gradient file":
    macro_gradMicrorotationFileName}}

    # Calls up the function to automatically select the boundary condi-
    # tions

    (bilinear_form, linear_form, boundary_conditions, 
    macro_quantitiesClasses, fields_namesDict, solution_fields, 
    variation_fields, trial_functions, monolithic_solution, 
    mixed_element, volume_inverse, fields_corrections) = multiscale_BCsTools.select_multiscaleBoundaryConditions(
    multiscale_BCsDict, elements_dictionary, mesh_dataClass, 
    fluctuation_field=False)

    ####################################################################
    #                         Variational forms                        #
    ####################################################################

    # Recovers the primal fields and their variations

    u_new = solution_fields["displacement"]
    
    phi_new = solution_fields["microrotation"]

    variation_u = variation_fields["displacement"]
    
    variation_phi = variation_fields["microrotation"]

    # Constructs the variational form for the inner work

    internal_VarForm = variational_tools.hyperelastic_micropolarInternalWorkFirstPiola(
    u_new, phi_new, variation_u, variation_phi, constitutive_model, 
    mesh_dataClass)

    # Constructs the residual and evaluates its derivative w.r.t. the
    # trial solution

    residual_form = ((volume_inverse*internal_VarForm)+bilinear_form-
    linear_form)

    residual_derivative = derivative(residual_form, monolithic_solution, 
    trial_functions)

    # Makes the boundary condition an empty list, because the macrosca-
    # le boundary conditions are applied directly onto the variational
    # form, using Lagrange multipliers

    Res = NonlinearVariationalProblem(residual_form, monolithic_solution,
    boundary_conditions, J=residual_derivative)

    ####################################################################
    #                    Solver parameters setting                     #
    ####################################################################

    solver = NonlinearVariationalSolver(Res)

    ####################################################################
    #                 Solution and pseudotime stepping                 #
    ####################################################################

    # Iterates through the pseudotime stepping algortihm 

    if len(solution_name)==0:

        for field_name in fields_namesDict:

            solution_name.append([field_name, "Microscale"])

    newton_raphson_tools.newton_raphsonMultipleFields(
    maximum_loadingSteps, solver, monolithic_solution, fields_namesDict, 
    mixed_element, mesh_dataClass, constitutive_model, 
    post_processesList=post_processes, solver_parameters=
    solver_parameters, solution_name=solution_name, 
    macro_quantitiesClasses=macro_quantitiesClasses, fields_corrections=
    fields_corrections)