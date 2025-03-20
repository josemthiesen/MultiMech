# Routine to store methods for pseudotime stepping

from dolfin import *

# Defines a function to iterate through a Newton-Raphson loop of a vari-
# ational problem of a single field

def newton_raphsonSingleField(t, t_final, delta_t, maximum_loadingSteps,
solver, solution_field, solution_fieldFile, dirichlet_loads=[],
neumann_loads=[], solver_parameters=dict()):
    
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

        print("\n#########################################################"+
        "###############\n#                 Incremental step: "+str(
        time_counter+1)+"; current time: "+str(t)+"               #\n#####"+
        "#################################################################"+
        "##\n")

        # Solves the nonlinear variational problem 

        solver.solve()

        solution_field.rename("DNS Displacement", "DNS")

        # Updates the files

        solution_fieldFile << solution_field

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
#                       Solver parameters
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

        solver.parameters["nonlinear_solver"] = "newton"

    else:

        solver.parameters["nonlinear_solver"] = solver_parameters["non"+
        "linear_solver"]

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
        'tolerance'] = parameter_types["krylov_relative_tolerance"]

    if "krylov_maximum_iterations" in parameter_types:

        solver.parameters['newton_solver']['krylov_solver']['maximum_i'+
        'terations'] = parameter_types["krylov_maximum_iterations"]

    if "krylov_monitor_convergence" in parameter_types:

        solver.parameters['newton_solver']['krylov_solver']['monitor_c'+
        'onvergence'] = parameter_types["krylov_monitor_convergence"]

    # Returns the updated solver

    return solver