# Routine to store methods for pseudotime stepping

from dolfin import *

import numpy as np

import source.tool_box.mesh_handling_tools as mesh_tools

########################################################################
#                        Newton-Raphson schemes                        #
########################################################################

# Defines a function to iterate through a Newton-Raphson loop of a vari-
# ational problem of a single field

def newton_raphsonSingleField(t, t_final, delta_t, maximum_loadingSteps,
solver, solution_field, post_processes, domain_meshCollection, 
post_processesSubmesh=dict(), dirichlet_loads=[], neumann_loads=[], 
solver_parameters=dict(), solution_name=["solution", "DNS"], 
volume_physGroupsSubmesh=[]):
    
    # Verifies if the physical groups for the submesh is an integer

    if isinstance(volume_physGroupsSubmesh, int):

        # Transforms into a list

        volume_physGroupsSubmesh = [volume_physGroupsSubmesh]
    
    # If there are volumetric physical groups to build a submesh

    if len(volume_physGroupsSubmesh)>0:

        (RVE_submesh, domain_meshFunction, UV_submesh, RVE_meshMapping, 
        parent_meshMapping, solution_submesh, RVE_toParentCellMap) = mesh_tools.create_submesh(
        mesh, domain_meshCollection, volume_physGroupsSubmesh, 
        solution_field)
    
    # Initializes a dictionary of post processes objects, files for e-
    # xample

    post_processingObjects = dict()

    for post_processName, post_process in post_processes.items():

        post_processingObjects[post_processName] = post_process[0](
        post_process[2], post_process[3])

    # Initializes a dictionary of post processes objects for the submesh, 
    # files for example, for each field

    post_processingObjectsSubmesh = dict()

    for post_processName, post_process in post_processesSubmesh.items():

        post_processingObjectsSubmesh[post_processName] = post_process[
        0](post_process[2], post_process[3])
    
    # Updates the solver parameters

    solver = set_solverParameters(solver, solver_parameters)
    
    # Verifies if there are no loads

    if len(dirichlet_loads)==0 and len(neumann_loads)==0:

        print("\nWARNING: there are no Dirichlet boundary conditions n"+
        "or Neumann boundary conditions\n")

    # Initializes the pseudotime counter

    time_counter = 0

    # Iterates through the pseudotime stepping

    while t<t_final:

        # Prints step information

        print_stepInfo(time_counter+1, t)

        # Solves the nonlinear variational problem 

        solver.solve()

        solution_field.rename(*solution_name)

        # Updates the post processes objects

        for post_processName, post_process in post_processes.items():

            post_processingObjects[post_processName] = post_process[1](
            post_processingObjects[post_processName], solution_field, t)

        # If a submesh is to be populated with part of the solution

        if len(volume_physGroupsSubmesh)>0 and (len(list(
        post_processesSubmesh.keys()))>0):

            solution_submesh = mesh_tools.field_parentToSubmesh(
            RVE_submesh, solution_submesh, solution_field, 
            RVE_toParentCellMap, RVE_meshMapping, parent_meshMapping)

            solution_submesh.rename(*solution_name)

            # Updates the post processes objects

            for post_processName, post_process in post_processesSubmesh.items():

                post_processingObjectsSubmesh[post_processName] = post_process[
                1](post_processingObjectsSubmesh[post_processName], 
                solution_submesh, t)

        # Updates the pseudo time variables and the counter

        t += delta_t
        
        time_counter += 1

        # Updates the Dirichlet boundary conditions 

        for dirichlet_load in dirichlet_loads:

            dirichlet_load.t = t

        # Updates the Neumann boundary conditions 

        for neumann_load in neumann_loads:

            neumann_load.t = t

        # Verifies if the maximum number of laoding steps has been 
        # reached

        if time_counter>=maximum_loadingSteps:

            print("\nThe maximum number of loading steps,",
            maximum_loadingSteps, "has just been reached. Stops the si"+
            "mulation immediatly\n")

            break

# Defines a function to iterate through a Newton-Raphson loop of a vari-
# ational problem of multiple fields

def newton_raphsonMultipleFields(t, t_final, delta_t, 
maximum_loadingSteps, solver, solution_field, post_processes, 
mixed_element, domain_meshCollection, post_processesSubmesh=[], 
dirichlet_loads=[], neumann_loads=[], solver_parameters=dict(), 
solution_name=[["solution", "DNS"]], volume_physGroupsSubmesh=[]):
    
    # Verifies if the physical groups for the submesh is an integer

    if isinstance(volume_physGroupsSubmesh, int):

        # Transforms into a list

        volume_physGroupsSubmesh = [volume_physGroupsSubmesh]
    
    # If there are volumetric physical groups to build a submesh

    if len(volume_physGroupsSubmesh)>0:

        (RVE_submesh, domain_meshFunction, UV_submesh, RVE_meshMapping, 
        parent_meshMapping, solution_submesh, RVE_toParentCellMap) = mesh_tools.create_submesh(
        mesh, domain_meshCollection, volume_physGroupsSubmesh, 
        solution_field, mixed_element=mixed_element)
    
    # Gets the number of fields in the mixed element

    n_fields = mixed_element.num_sub_elements()

    # Verifies if the post processes is a list

    if not isinstance(post_processes, list):

        raise ValueError("post_processes must be a list of dictionarie"+
        "s in newton_raphsonMultipleFields, because there must be post"+
        "-processing steps for each field.")
    
    # Initializes a dictionary of post processes objects, files for e-
    # xample, for each field

    post_processingObjects = []

    for i in range(n_fields):

        post_processingObjects.append(dict())

        for post_processName, post_process in post_processes[i].items():

            post_processingObjects[-1][post_processName] = post_process[
            0](post_process[2], post_process[3])

    # Makes the same thing for the submesh post-processes

    if not isinstance(post_processesSubmesh, list):

        raise ValueError("post_processesSubmesh must be a list of dict"+
        "ionaries in newton_raphsonMultipleFields, because there must "+
        "be post-processing steps for each field.")
    
    # Initializes a dictionary of post processes objects, files for e-
    # xample, for each field

    post_processingObjectsSubmesh = []

    for i in range(n_fields):

        post_processingObjectsSubmesh.append(dict())

        for post_processName, post_process in post_processesSubmesh[i
        ].items():

            post_processingObjectsSubmesh[-1][post_processName] = post_process[
            0](post_process[2], post_process[3])
    
    # Updates the solver parameters

    solver = set_solverParameters(solver, solver_parameters)
    
    # Verifies if there are no loads

    if len(dirichlet_loads)==0 and len(neumann_loads)==0:

        print("\nWARNING: there are no Dirichlet boundary conditions n"+
        "or Neumann boundary conditions\n")

    # Initializes the pseudotime counter

    time_counter = 0

    # Iterates through the pseudotime stepping

    while t<t_final:

        # Prints step information

        print_stepInfo(time_counter+1, t)

        # Solves the nonlinear variational problem 

        solver.solve()

        # Splits the solution

        split_solution = list(solution_field.split(deepcopy=True))

        # Post-process the solution

        if len(post_processes)>0:

            for i in range(n_fields):

                # Updates the name

                if len(solution_name)==n_fields:

                    split_solution[i].rename(*solution_name[i])

                # Updates the post processes objects

                for post_processName, post_process in post_processes[i
                ].items():

                    post_processingObjects[i][post_processName] = post_process[
                    1](post_processingObjects[i][post_processName], 
                    split_solution[i], t)

        # If a submesh is to be populated with part of the solution

        if len(volume_physGroupsSubmesh)>0 and (len(post_processesSubmesh
        )>0):

            solution_submesh = mesh_tools.field_parentToSubmesh(
            RVE_submesh, solution_submesh, solution_field, 
            RVE_toParentCellMap, RVE_meshMapping, parent_meshMapping)

            # Splits the solution

            split_solutionSubmesh = list(solution_submesh.split(deepcopy
            =True))

            # Post-process the submesh solution

            for i in range(n_fields):

                # Updates the name of the fields

                if len(solution_name)==n_fields:

                    split_solutionSubmesh[i].rename(*solution_name[i])

                # Updates the post processes objects

                for post_processName, post_process in post_processesSubmesh[
                i].items():

                    post_processingObjectsSubmesh[i][post_processName] = post_process[
                    1](post_processingObjectsSubmesh[i][post_processName], 
                    split_solutionSubmesh[i], t)

        # Updates the pseudo time variables and the counter

        t += delta_t
        
        time_counter += 1

        # Updates the Dirichlet boundary conditions 

        for dirichlet_load in dirichlet_loads:

            dirichlet_load.t = t

        # Updates the Neumann boundary conditions 

        for neumann_load in neumann_loads:

            neumann_load.t = t

        # Verifies if the maximum number of laoding steps has been 
        # reached

        if time_counter>=maximum_loadingSteps:

            print("\nThe maximum number of loading steps,",
            maximum_loadingSteps, "has just been reached. Stops the si"+
            "mulation immediatly\n")

            break

########################################################################
#                              Utilities                               #
########################################################################

# Defines a function to update solver parameters

def set_solverParameters(solver, solver_parameters):

    # Sets a list of implemented solver parameters

    admissible_keys = ["nonlinear_solver", "linear_solver", "newton_re"+
    "lative_tolerance", "newton_absolute_tolerance", "newton_maximum_i"+
    "terations", "preconditioner", "krylov_absolute_tolerance", "krylo"+
    "v_relative_tolerance", "krylov_maximum_iterations", "krylov_monit"+
    "or_convergence"]

    # Gets the keys of the solver parameters dictionary

    parameter_types = solver_parameters.keys()

    # Iterates the keys of the solver parameters to verify if any of 
    # them is not admissible

    for key in parameter_types:

        if not (key in admissible_keys):

            raise NameError("The key "+str(key)+" is not an admissible"+
            " key to set solver parameters.")
        
    # Sets the solver parameters

    if "nonlinear_solver" in parameter_types:

        solver.parameters["nonlinear_solver"] = solver_parameters["non"+
        "linear_solver"]

    else:

        solver.parameters["nonlinear_solver"] = "newton"

    if "linear_solver" in parameter_types:

        solver.parameters["newton_solver"]["linear_solver"] = (
        solver_parameters["linear_solver"])

    if "newton_relative_tolerance" in parameter_types:

        solver.parameters["newton_solver"]["relative_tolerance"] = (
        solver_parameters["newton_relative_tolerance"])

    if "newton_absolute_tolerance" in parameter_types:

        solver.parameters["newton_solver"]["absolute_tolerance"] = (
        solver_parameters["newton_absolute_tolerance"])

    if "newton_maximum_iterations" in parameter_types:

        solver.parameters["newton_solver"]["maximum_iterations"] = (
        solver_parameters["newton_maximum_iterations"])

    if "preconditioner" in parameter_types:

        solver.parameters["newton_solver"]["preconditioner"] = (
        solver_parameters["preconditioner"])

    if "krylov_absolute_tolerance" in parameter_types:

        solver.parameters['newton_solver']['krylov_solver']['absolute_'+
        'tolerance'] = solver_parameters["krylov_absolute_tolerance"]

    if "krylov_relative_tolerance" in parameter_types:

        solver.parameters['newton_solver']['krylov_solver']['relative_'+
        'tolerance'] = solver_parameters["krylov_relative_tolerance"]

    if "krylov_maximum_iterations" in parameter_types:

        solver.parameters['newton_solver']['krylov_solver']['maximum_i'+
        'terations'] = solver_parameters["krylov_maximum_iterations"]

    if "krylov_monitor_convergence" in parameter_types:

        solver.parameters['newton_solver']['krylov_solver']['monitor_c'+
        'onvergence'] = solver_parameters["krylov_monitor_convergence"]

    # Returns the updated solver

    return solver

# Defines a function to print the stepping information

def print_stepInfo(step, time):

    # Transforms the time and the step into strings

    step = str(step)

    time = str(time)

    # Evaluates the number of characters 

    n_characters = len(step)+len(time)

    # Verifies if there is enough space for the whole step information

    if n_characters<=36:

        # Calculates the clearance to each side

        clearance_left = int(np.ceil(0.5*(36-n_characters)))

        clearance_right = 36-clearance_left-n_characters

        # Makes the clearance space

        clearance_spaceLeft = ""

        clearance_spaceRight = ""

        for i in range(clearance_left):

            clearance_spaceLeft += " "

        for i in range(clearance_right):

            clearance_spaceRight += " "

        print("\n#####################################################"+
        "###################\n#"+clearance_spaceLeft+"Incremental step"+
        ": "+step+"; current time: "+time+clearance_spaceRight+"#\n###"+
        "#############################################################"+
        "########\n")

    else:

        # Calculates the clearance to each side

        clearance_left = int(np.ceil(0.5*(52-len(step))))

        clearance_right = 52-clearance_left-len(step)

        # Makes the clearance space

        clearance_spaceLeft = ""

        clearance_spaceRight = ""

        for i in range(clearance_left):

            clearance_spaceLeft += " "

        for i in range(clearance_right):

            clearance_spaceRight += " "

        print("\n#####################################################"+
        "###################\n#"+clearance_spaceLeft+"Incremental step"+
        ": "+step+clearance_spaceRight+"#")

        # Calculates the clearance to each side

        clearance_left = int(np.ceil(0.5*(56-len(time))))

        clearance_right = 56-clearance_left-len(time)

        # Makes the clearance space

        clearance_spaceLeft = ""

        clearance_spaceRight = ""

        for i in range(clearance_left):

            clearance_spaceLeft += " "

        for i in range(clearance_right):

            clearance_spaceRight += " "

        print("#"+clearance_spaceLeft+"current time: "+time+
        clearance_spaceRight+"#\n#####################################"+
        "###################################\n")