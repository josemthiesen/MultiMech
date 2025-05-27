# Routine to store the variational form and other accessories for a hy-
# perelastic micropolar (Cosserat) continuum in solid mechanics

from dolfin import *

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.functional_tools as functional_tools

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

    # Assembles a dictionary of finite elements for the two primal 
    # fields: displacement and microrotation. Each field has a key and
    # the corresponding value is another dictionary, which has keys for
    # necessary information to create finite elements

    elements_dictionary = {"displacement": {"field type": "vector", "i"+
    "nterpolation function": "CG", "polynomial degree": 
    polynomial_degreeDisplacement}, "microrotation": {"field type": "v"+
    "ector", "interpolation function": "CG", "polynomial degree": 
    polynomial_degreeMicrorotation}}

    # From the dictionary of elements, the finite elements are created,
    # then, all the rest is created: the function spaces, trial and test 
    # functions, solution function. Everything is split and named by ac-
    # cording to the element's name

    (monolithic_functionSpace, solution_new, fields_names, 
    solution_fields, variation_fields, delta_solution, fields_namesDict
    ) = functional_tools.construct_monolithicFunctionSpace(
    elements_dictionary, mesh_dataClass, verbose=verbose)

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

    # Gets the individual solution fields and the individual variations
    # of the fields

    u_new = solution_fields["displacement"]

    phi_new = solution_fields["microrotation"]

    variation_u = variation_fields["displacement"]
    
    variation_phi = variation_fields["microrotation"]

    # Constructs the variational form for the inner work

    internal_VarForm = variational_tools.hyperelastic_micropolarInternalWorkFirstPiola(
    u_new, phi_new, variation_u, variation_phi, constitutive_model, 
    mesh_dataClass)

    # Constructs the variational forms for the traction work

    traction_VarForm, neumann_loads = variational_tools.traction_work(
    traction_dictionary, u_new, variation_u, mesh_dataClass, 
    neumann_loads)

    # Constructs the variational forms for the moment work on the boun-
    # dary. Note that the function traction_work was reused, because the
    # variational construction is the same for traction and for moment

    moment_VarForm, neumann_loads = variational_tools.traction_work(
    moment_dictionary, phi_new, variation_phi, mesh_dataClass, 
    neumann_loads)

    ####################################################################
    #              Problem and solver parameters setting               #
    ####################################################################

    # Assembles the residual and the nonlinear problem object. Sets the
    # solver parameters too

    residual_form = internal_VarForm-traction_VarForm-moment_VarForm

    solver = functional_tools.set_nonlinearProblem(residual_form, 
    solution_new, delta_solution, bc, solver_parameters=
    solver_parameters)

    ####################################################################
    #                 Solution and pseudotime stepping                 #
    ####################################################################

    # Iterates through the pseudotime stepping algorithm 

    if len(solution_name)==0:

        for field_name in fields_namesDict:

            solution_name.append([field_name, "DNS"])

    newton_raphson_tools.newton_raphsonMultipleFields(solver, 
    solution_new, fields_namesDict, mesh_dataClass, constitutive_model, 
    post_processesList=post_processes, post_processesSubmeshList=
    post_processesSubmesh, dirichlet_loads=dirichlet_loads, neumann_loads
    =neumann_loads, volume_physGroupsSubmesh=volume_physGroupsSubmesh, 
    solution_name=solution_name, t=t, t_final=t_final, 
    maximum_loadingSteps=maximum_loadingSteps)