# Routine to store the variational form and other accessories for a hy-
# perelastic micropolar (Cosserat) continuum in solid mechanics

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.pseudotime_stepping_tools as newton_raphson_tools

import source.tool_box.programming_tools as programming_tools

# Defines a function to model a hyperelastic problem with a displacement
# and a microrotation fields only

@programming_tools.optional_argumentsInitializer({'neumann_loads': 
lambda: [], 'dirichlet_loads': lambda: [], 'solution_name': lambda: [], 
'simple_supportDisplacementPhysicalGroups': lambda: dict(), ('simple_s'+
'upportMicrorotationPhysicalGroups'): lambda: dict(), ('volume_physGro'+
'upsSubmesh'): lambda: [], 'post_processesSubmesh': lambda: []})

def hyperelasticity_displacementMicrorotationBased(constitutive_model, 
traction_dictionary, moment_dictionary, maximum_loadingSteps, t_final, 
post_processes, mesh_fileName, solver_parameters, 
polynomial_degreeDisplacement=2, polynomial_degreeMicrorotation=2, 
t=0.0, fixed_supportDisplacementPhysicalGroups=0, neumann_loads=None, 
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

    mixed_element = MixedElement([displacement_element, 
    microrotation_element])

    # Defines the finite element space for the monolithic solution

    monolithic_functionSpace = FunctionSpace(mesh_dataClass.mesh, 
    mixed_element)

    # Sets the names for each field that will be used to retrive and 
    # sort post-processes into a dictionary. The names are the keys and
    # the values are the indexes in the monolithic function space

    fields_names = {"displacement":0, "microrotation":1}

    if verbose:

        print("Finishes creating the function spaces\n")

    ####################################################################
    #                        Boundary conditions                       #
    ####################################################################

    # Defines the boundary conditions for fixed facets in terms of dis-
    # placement

    bc = BCs_tools.fixed_supportDirichletBC(monolithic_functionSpace, 
    mesh_dataClass, boundary_physicalGroups=
    fixed_supportDisplacementPhysicalGroups, sub_fieldsToApplyBC=[0])

    # Defines the boundary conditions for fixed facets in terms of mi-
    # crorotation

    bc = BCs_tools.fixed_supportDirichletBC(monolithic_functionSpace, 
    mesh_dataClass, boundary_physicalGroups=
    fixed_supportMicrorotationPhysicalGroups, sub_fieldsToApplyBC=[1],
    boundary_conditions=bc)

    # Adds boundary conditions for simply supported facets in terms of
    # displacement

    bc = BCs_tools.simple_supportDirichletBC(monolithic_functionSpace, 
    mesh_dataClass, boundary_physicalGroups=
    simple_supportDisplacementPhysicalGroups, boundary_conditions=bc,
    sub_fieldsToApplyBC=[0])

    # Adds boundary conditions for simply supported facets in terms of
    # microrotation

    bc = BCs_tools.simple_supportDirichletBC(monolithic_functionSpace, 
    mesh_dataClass, boundary_physicalGroups=
    simple_supportMicrorotationPhysicalGroups, boundary_conditions=bc,
    sub_fieldsToApplyBC=[1])

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

    #traction_VarForm = (dot(as_vector([0.0, neumann_loads[0], 0.0]), 
    #variation_u)*ds(6))

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

    Res = NonlinearVariationalProblem(residual_form, solution_new, bc, 
    J=residual_derivative)

    ####################################################################
    #                    Solver parameters setting                     #
    ####################################################################

    solver = NonlinearVariationalSolver(Res)

    ####################################################################
    #                 Solution and pseudotime stepping                 #
    ####################################################################

    # Iterates through the pseudotime stepping algorithm 

    if len(solution_name)==0:

        for field_name in fields_names:

            solution_name.append([field_name, "DNS"])

    newton_raphson_tools.newton_raphsonMultipleFields(
    maximum_loadingSteps, solver, solution_new, fields_names,
    mixed_element, mesh_dataClass, constitutive_model, 
    post_processesList=post_processes, post_processesSubmeshList=
    post_processesSubmesh, dirichlet_loads=dirichlet_loads, 
    neumann_loads=neumann_loads, solver_parameters=solver_parameters, 
    volume_physGroupsSubmesh=volume_physGroupsSubmesh, solution_name=
    solution_name, t=t, t_final=t_final)