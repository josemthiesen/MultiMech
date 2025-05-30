# Routine to store the variational form and other accessories for a hy-
# perelastic Cauchy continuum in solid mechanics

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.functional_tools as functional_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

import source.tool_box.programming_tools as programming_tools

# Defines a function to model a hyperelastic problem with a displacement
# field only

@programming_tools.optional_argumentsInitializer({'neumann_loads': 
lambda: [], 'dirichlet_loads': lambda: [], 'solution_name': lambda: [
"solution", "DNS"], 'simple_supportPhysicalGroups': lambda: dict(), 
'volume_physGroupsSubmesh': lambda: [], 'post_processesSubmesh': lambda: 
dict(), 'prescribed_displacement': lambda: dict()})

def hyperelasticity_displacementBased(constitutive_model, 
traction_dictionary, maximum_loadingSteps, t_final, post_processes, 
mesh_fileName, solver_parameters, neumann_loads=None, dirichlet_loads=
None, prescribed_displacement=None, polynomial_degree=2, 
quadrature_degree=2, t=0.0, fixed_supportPhysicalGroups=0, 
simple_supportPhysicalGroups=None, volume_physGroupsSubmesh=None, 
post_processesSubmesh=None, solution_name=None, verbose=False):

    ####################################################################
    #                               Mesh                               #
    ####################################################################

    # Reads the mesh and constructs some fenics objects using the xdmf 
    # file

    mesh_dataClass = mesh_tools.read_mshMesh(mesh_fileName, 
    quadrature_degree=quadrature_degree, verbose=verbose)

    ####################################################################
    #                          Function space                          #
    ####################################################################

    # Assembles a dictionary of finite elements for the two primal 
    # fields: displacement and microrotation. Each field has a key and
    # the corresponding value is another dictionary, which has keys for
    # necessary information to create finite elements

    elements_dictionary = {"displacement": {"field type": "vector", "i"+
    "nterpolation function": "CG", "polynomial degree": 
    polynomial_degree}}

    # From the dictionary of elements, the finite elements are created,
    # then, all the rest is created: the function spaces, trial and test 
    # functions, solution function. Everything is split and named by ac-
    # cording to the element's name

    (solution_functionSpace, solution_new, fields_names, solution_fields, 
    variation_fields, delta_solution, fields_namesDict
    ) = functional_tools.construct_monolithicFunctionSpace(
    elements_dictionary, mesh_dataClass, verbose=verbose)

    ####################################################################
    #                        Boundary conditions                       #
    ####################################################################

    # Defines the boundary conditions for fixed facets

    bc = BCs_tools.fixed_supportDirichletBC(solution_functionSpace, 
    mesh_dataClass, boundary_physicalGroups=fixed_supportPhysicalGroups,
    boundary_conditions=[])

    # Adds boundary conditions for simply supported facets

    bc = BCs_tools.simple_supportDirichletBC(solution_functionSpace, 
    mesh_dataClass, simple_supportPhysicalGroups, boundary_conditions=bc)

    # Adds prescribed displacement using the Dirichlet loads

    bc = BCs_tools.prescribed_DirichletBC(prescribed_displacement, 
    solution_functionSpace, mesh_dataClass, boundary_conditions=bc)

    ####################################################################
    #                         Variational forms                        #
    ####################################################################

    # Constructs the variational form for the inner work

    internal_VarForm = variational_tools.hyperelastic_internalWorkFirstPiola(
    "displacement", solution_fields, variation_fields, 
    constitutive_model, mesh_dataClass)

    # Constructs the variational forms for the traction work

    traction_VarForm, neumann_loads = variational_tools.traction_work(
    traction_dictionary, "displacement", solution_fields, 
    variation_fields, solution_new, fields_namesDict, mesh_dataClass, 
    neumann_loads)

    ####################################################################
    #              Problem and solver parameters setting               #
    ####################################################################

    # Assembles the residual and the nonlinear problem object. Sets the
    # solver parameters too

    residual_form = internal_VarForm-traction_VarForm

    solver = functional_tools.set_nonlinearProblem(residual_form, 
    solution_new, delta_solution, bc, solver_parameters=
    solver_parameters)

    ####################################################################
    #                 Solution and pseudotime stepping                 #
    ####################################################################

    # Iterates through the pseudotime stepping algortihm 

    newton_raphson_tools.newton_raphsonSingleField(solver, solution_new, 
    fields_namesDict, mesh_dataClass, constitutive_model, 
    post_processesDict=post_processes, post_processesSubmeshDict=
    post_processesSubmesh, neumann_loads=neumann_loads, dirichlet_loads=
    dirichlet_loads, solution_name=solution_name, 
    volume_physGroupsSubmesh=volume_physGroupsSubmesh, t=t, t_final=
    t_final, maximum_loadingSteps=maximum_loadingSteps)