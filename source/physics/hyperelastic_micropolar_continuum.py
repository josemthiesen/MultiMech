# Routine to store the variational form and other accessories for a hy-
# perelastic micropolar (Cosserat) continuum in solid mechanics

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

# Defines a function to model a hyperelastic problem with a displacement
# and a microrotation fields only

def hyperelasticity_displacementMicrorotationBased(constitutive_model, 
traction_dictionary, neumann_loads, maximum_loadingSteps, t_final, 
post_processes, mesh_fileName, solver_parameters, polynomial_degree=2, 
t=0.0, fixed_supportPhysicalGroups=0, simple_supportPhysicalGroups=dict(
), volume_physGroupsSubmesh=[], post_processesSubmesh=dict()):

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

    # Constructs elements for the displacement and for the microrotation
    # fields

    displacement_element = VectorElement("CG", mesh.ufl_cell(), 
    polynomial_degree)

    microrotation_element = VectorElement("CG", mesh.ufl_cell(), 
    polynomial_degree)

    mixed_element = MixedElement([displacement_element, 
    microrotation_element])

    # Defines the finite element space for the monolithic solution

    monolithic_functionSpace = FunctionSpace(mesh, mixed_element)

    ####################################################################
    #                        Boundary conditions                       #
    ####################################################################

    # Defines the boundary conditions for fixed facets

    bc = BCs_tools.fixed_supportDirichletBC(monolithic_functionSpace, 
    boundary_meshFunction, boundary_physicalGroups=
    fixed_supportPhysicalGroups)

    # Adds boundary conditions for simply supported facets

    bc = BCs_tools.simple_supportDirichletBC(monolithic_functionSpace, 
    boundary_meshFunction, boundary_physicalGroups=
    simple_supportPhysicalGroups, boundary_conditions=bc)

    ####################################################################
    #                         Variational forms                        #
    ####################################################################

    # Defines the trial and test functions

    delta_solution = TrialFunction(monolithic_functionSpace) 

    variation_solution = TestFunction(monolithic_functionSpace)

    # Creates the function for the updated solution, i.e. the vector of 
    # parameters

    solution_new = Function(monolithic_functionSpace)

    # Constructs the variational form for the inner work

    internal_VarForm = variational_tools.hyperelastic_internalWorkFirstPiola(
    solution_new, variation_solution, constitutive_model, dx)

    # Constructs the variational forms for the traction work

    traction_VarForm = variational_tools.traction_work(
    traction_dictionary, variation_solution, ds)

    # Assembles the residual, takes the Gateaux derivative and assembles
    # the nonlinear problem object

    residual_form = internal_VarForm-traction_VarForm

    residual_derivative = derivative(residual_form, solution_new, 
    delta_solution)

    Res = NonlinearVariationalProblem(residual_form, solution_new, bc, 
    J=residual_derivative)

    ####################################################################
    #                    Solver parameters setting                     #
    ####################################################################

    solver = NonlinearVariationalSolver(Res)

    ####################################################################
    #                 Solution and pseudotime stepping                 #
    ####################################################################

    # Evaluates the pseudotime step

    delta_t = (t_final-t)/maximum_loadingSteps

    # Iterates through the pseudotime stepping algortihm 

    newton_raphson_tools.newton_raphsonSingleField(t, t_final, delta_t, 
    maximum_loadingSteps, solver, solution_new, domain_meshCollection, 
    constitutive_model, dx, post_processesDict=post_processes, 
    post_processesSubmeshDict=post_processesSubmesh, neumann_loads=
    neumann_loads, solver_parameters=solver_parameters, 
    volume_physGroupsSubmesh=volume_physGroupsSubmesh)