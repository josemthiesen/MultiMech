# Routine to store the variational form and other accessories for a hy-
# perelastic solid mechanics problem

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.file_handling_tools as file_handling_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

# Defines a function to model a hyperelastic problem with a displacement
# field only

def hyperelasticity_displacementBased(constitutive_model, 
traction_dictionary, neumann_loads, maximum_loadingSteps, t_final, 
results_path, displacement_fileName, mesh_fileName, solver_parameters, 
polynomial_degree=2, t= 0.0, fixed_supportPhysicalGroups=0, 
simple_supportPhysicalGroups=dict()):

    ####################################################################
    #                               Mesh                               #
    ####################################################################

    # Reads the mesh and constructs some fenics objects using the xdmf 
    # file

    (mesh, dx, ds, n, domain_meshCollection, domain_meshFunction, 
    boundary_meshCollection, boundary_meshFunction) = mesh_tools.read_mshMesh(
    mesh_fileName)

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
    fixed_supportPhysicalGroups)

    # Adds boundary conditions for simply supported facets

    bc = BCs_tools.simple_supportDirichletBC(U, boundary_meshFunction,
    simple_supportPhysicalGroups, boundary_conditions=bc)

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
    u_new, v, constitutive_model, dx)

    # Constructs the variational forms for the traction work

    traction_VarForm = variational_tools.traction_work(
    traction_dictionary, v, ds)

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
    #                       Files initialization                       #
    ####################################################################

    # Creates the path to the displacement file

    displacement_file = file_handling_tools.verify_path(results_path, 
    displacement_fileName)

    displacement_file = File(displacement_file)

    ####################################################################
    #                 Solution and pseudotime stepping                 #
    ####################################################################

    # Evaluates the pseudotime step

    delta_t = (t_final-t)/maximum_loadingSteps

    # Iterates through the pseudotime stepping algortihm 

    newton_raphson_tools.newton_raphsonSingleField(t, t_final, delta_t, 
    maximum_loadingSteps, solver, u_new, displacement_file,neumann_loads=
    neumann_loads, solver_parameters=solver_parameters)