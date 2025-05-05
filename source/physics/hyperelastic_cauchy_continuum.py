# Routine to store the variational form and other accessories for a hy-
# perelastic Cauchy continuum in solid mechanics

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

import source.tool_box.variational_tools as variational_tools

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

    # Defines the finite element spaces for the displacement field, u

    print("Polynomial degree:", polynomial_degree, "\n")

    U = VectorFunctionSpace(mesh_dataClass.mesh, "Lagrange", 
    polynomial_degree)

    ####################################################################
    #                        Boundary conditions                       #
    ####################################################################

    # Defines the boundary conditions for fixed facets

    bc = BCs_tools.fixed_supportDirichletBC(U, mesh_dataClass, 
    boundary_physicalGroups=fixed_supportPhysicalGroups,
    boundary_conditions=[])

    # Adds boundary conditions for simply supported facets

    bc = BCs_tools.simple_supportDirichletBC(U, mesh_dataClass,
    simple_supportPhysicalGroups, boundary_conditions=bc)

    # Adds prescribed displacement using the Dirichlet loads

    bc = BCs_tools.prescribed_DirichletBC(prescribed_displacement, U,
    mesh_dataClass, boundary_conditions=bc)

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
    u_new, v, constitutive_model, mesh_dataClass)

    #Piola = constitutive_model.first_piolaStress(u_new)

    #internal_VarForm = (inner(Piola, grad(v))*dx)

    # Constructs the variational forms for the traction work

    traction_VarForm = variational_tools.traction_work(
    traction_dictionary, v, mesh_dataClass)

    #traction_VarForm = dot(as_vector([0.0, neumann_loads[0], 0.0]), v)*ds(6)

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
    maximum_loadingSteps, solver, u_new, mesh_dataClass, 
    constitutive_model, post_processesDict=post_processes, 
    post_processesSubmeshDict=post_processesSubmesh, neumann_loads=
    neumann_loads, dirichlet_loads=dirichlet_loads, solver_parameters=
    solver_parameters, volume_physGroupsSubmesh=volume_physGroupsSubmesh,
    solution_name=solution_name)