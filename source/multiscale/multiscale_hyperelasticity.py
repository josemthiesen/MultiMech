# Routine to store the multiscale methods for micropolar physics

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

import source.tool_box.programming_tools as programming_tools

import source.tool_box.multiscale_boundary_conditions_tools as multiscale_BCsTools

# Defines a function to model a hyperelastic problem with a displacement
# only in the microscale. It uses the macro quantities read from txt fi-
# les

@programming_tools.optional_argumentsInitializer({'neumann_loads': 
lambda: [], 'dirichlet_loads': lambda: [], 'solution_name': lambda: []})

def hyperelastic_microscale(displacement_multiscaleBC, 
macro_displacementFileName, macro_gradDisplacementFileName, 
constitutive_model, maximum_loadingSteps, post_processes, mesh_fileName, 
solver_parameters, polynomial_degree=2, solution_name=None, verbose=
False, fluctuation_field=False):

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

    # Assembles a dictionary of finite elements for the one and primal 
    # field: displacement. The field has a key and the corresponding va-
    # lue is another dictionary, which has keys for necessary informa-
    # tion to create finite elements

    elements_dictionary = {"displacement": {"field type": "vector", "i"+
    "nterpolation function": "CG", "polynomial degree": 
    polynomial_degree}}

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

    # Calls up the function to automatically select the boundary condi-
    # tions

    (bilinear_form, linear_form, boundary_conditions, 
    macro_quantitiesClasses, fields_namesDict, solution_fields, 
    variation_fields, trial_functions, monolithic_solution, 
    mixed_element, volume_inverse, fields_corrections) = multiscale_BCsTools.select_multiscaleBoundaryConditions(
    multiscale_BCsDict, elements_dictionary, mesh_dataClass, 
    fluctuation_field=fluctuation_field)

    ####################################################################
    #                         Variational forms                        #
    ####################################################################

    # Recovers the primal fields and their variations

    u_new = solution_fields["displacement"]

    variation_u = variation_fields["displacement"]

    # Constructs the variational form for the inner work

    internal_VarForm = variational_tools.hyperelastic_internalWorkFirstPiola(
    u_new, variation_u, constitutive_model, mesh_dataClass)

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

    # The hyperelastic microscale problem is, in nature, a single-field
    # problem. But, if minimally constrained boundary condition is used,
    # the problems gains one field for each Lagrange multiplier. Thus,
    # the solution of the problem is given to an algorithms that can 
    # handle multiple fields. 
    # If the solution, however, does not have multiple fields, the fol-
    # lowing function will automatically reconnect with the single-field
    # framework

    newton_raphson_tools.newton_raphsonMultipleFields(
    maximum_loadingSteps, solver, monolithic_solution, fields_namesDict, 
    mixed_element, mesh_dataClass, constitutive_model, 
    post_processesList=post_processes, solver_parameters=
    solver_parameters, solution_name=solution_name, 
    macro_quantitiesClasses=macro_quantitiesClasses, fields_corrections=
    fields_corrections)