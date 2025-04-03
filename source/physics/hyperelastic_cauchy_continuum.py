# Routine to store the variational form and other accessories for a hy-
# perelastic Cauchy continuum in solid mechanics

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

# Defines a function to model a hyperelastic problem with a displacement
# field only

def hyperelasticity_displacementBased(constitutive_model, 
traction_dictionary, maximum_loadingSteps, t_final, post_processes, 
mesh_fileName, solver_parameters, neumann_loads=[], dirichlet_loads=[],  
polynomial_degree=2, t=0.0, fixed_supportPhysicalGroups=0, 
simple_supportPhysicalGroups=dict(), volume_physGroupsSubmesh=[], 
post_processesSubmesh=dict(), solution_name=["solution", "DNS"]):

    ####################################################################
    #                               Mesh                               #
    ####################################################################

    # Reads the mesh and constructs some fenics objects using the xdmf 
    # file

    (mesh, dx, ds, n, domain_meshCollection, domain_meshFunction, 
    boundary_meshCollection, boundary_meshFunction, 
    domain_physGroupsNamesToTags, boundary_physGroupsNamesToTags
    ) = mesh_tools.read_mshMesh(mesh_fileName)

    ####################################################################
    #                          Function space                          #
    ####################################################################

    # Defines the finite element spaces for the displacement field, u

    U = VectorFunctionSpace(mesh, "Lagrange", polynomial_degree)

    ####################################################################
    #                        Boundary conditions                       #
    ####################################################################

    # Defines the boundary conditions for fixed facets

    bc = BCs_tools.fixed_supportDirichletBC(U, boundary_meshFunction, 
    boundary_physicalGroups=fixed_supportPhysicalGroups,
    boundary_physGroupsNamesToTags=boundary_physGroupsNamesToTags)

    # Adds boundary conditions for simply supported facets

    bc = BCs_tools.simple_supportDirichletBC(U, boundary_meshFunction,
    simple_supportPhysicalGroups, boundary_conditions=bc,
    boundary_physGroupsNamesToTags=boundary_physGroupsNamesToTags)

    ####################################################################
    #                         Variational forms                        #
    ####################################################################

    # Defines the trial and test functions

    delta_u = TrialFunction(U) 

    v = TestFunction(U)

    # Creates the function for the updated solution, i.e. the vector of 
    # parameters

    u_new = Function(U)

    # Constructs the variational form for the inner work

    internal_VarForm = variational_tools.hyperelastic_internalWorkFirstPiola(
    u_new, v, constitutive_model, dx, domain_physGroupsNamesToTags=
    domain_physGroupsNamesToTags)

    # Constructs the variational forms for the traction work

    traction_VarForm = variational_tools.traction_work(
    traction_dictionary, v, ds, boundary_physGroupsNamesToTags=
    boundary_physGroupsNamesToTags)

    # Assembles the residual, takes the Gateaux derivative and assembles
    # the nonlinear problem object

    residual_form = internal_VarForm-traction_VarForm

    residual_derivative = derivative(residual_form , u_new, delta_u)

    Res = NonlinearVariationalProblem(residual_form, u_new, bc, J=
    residual_derivative)

    ####################################################################
    #                    Solver parameters setting                     #
    ####################################################################

    solver = NonlinearVariationalSolver(Res)

    ####################################################################
    #                 Solution and pseudotime stepping                 #
    ####################################################################

    # Iterates through the pseudotime stepping algortihm 

    newton_raphson_tools.newton_raphsonSingleField(t, t_final, 
    maximum_loadingSteps, solver, u_new, domain_meshCollection, 
    constitutive_model, dx, post_processesDict=post_processes, 
    post_processesSubmeshDict=post_processesSubmesh, neumann_loads=
    neumann_loads, dirichlet_loads=dirichlet_loads, solver_parameters=
    solver_parameters, volume_physGroupsSubmesh=volume_physGroupsSubmesh,
    solution_name=solution_name)