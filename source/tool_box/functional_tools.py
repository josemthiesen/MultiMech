# Routine to store functions to help with functional analysis tasks and
# programmatic tools

from dolfin import *

from abc import ABC, abstractmethod

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.numerical_tools as numerical_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.programming_tools as programming_tools

########################################################################
#                            Solver setting                            #
########################################################################

# Defines a function to set the nonlinear problem from the residual and
# the Gateaux derivative

def set_nonlinearProblem(residual_form, monolithic_solution, 
trial_functions, boundary_conditions, solver_parameters=None):
    
    # Evaluates the Gateaux derivative of the residue form
    
    residual_derivative = derivative(residual_form, monolithic_solution, 
    trial_functions)

    # Sets the nonlinear variational problem

    Res = NonlinearVariationalProblem(residual_form, monolithic_solution, 
    boundary_conditions, J=residual_derivative)

    # Sets the solver to this problem

    solver = NonlinearVariationalSolver(Res)

    if not (solver_parameters is None):

        solver = set_solverParameters(solver, solver_parameters)

    return solver

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

########################################################################
#                       Finite element creation                        #
########################################################################

# Defines a function to select a field from a list of fields

def select_fields(split_fieldsList, required_fieldsNames, 
fields_namesDict):
    
    # Tests the solution is a list

    if not isinstance(split_fieldsList, list):

        # If not, simply gives the field

        return split_fieldsList

    # Initializes a list of fields that will be used

    retrieved_fields = []

    # Iterates through the required fields' names

    for name in required_fieldsNames:

        # Verifies if this name is in the dictionary of fields' names 
        # that came from the generation of the finite elements

        if not (name in fields_namesDict):
                
            raise NameError("The field '"+str(name)+"' is required in "+
            "this post-process, but it is not in the list of names of "+
            "fields, that is: "+str(fields_namesDict.keys()))
        
        # Appends the field by its name

        retrieved_fields.append(split_fieldsList[fields_namesDict[name]])

    # Returns the list of fields

    if len(retrieved_fields)==1:

        return retrieved_fields[0]

    return retrieved_fields

# Defines a function to construct a function space given a dictionary of
# instructions to create finite elements

def construct_monolithicFunctionSpace(elements_dictionary, 
mesh_dataClass, verbose=False):
    
    # Transforms the dictionary of instructions into real finite elements
    # and get the names of the fields

    elements_dictionary, fields_names = construct_elementsDictionary(
    elements_dictionary, mesh_dataClass)

    # Constructs the mixed element

    mixed_element = 0

    if len(fields_names)>1:

        mixed_element = MixedElement([elements_dictionary[field_name
        ] for field_name in (fields_names)])

    else:

        mixed_element = elements_dictionary[fields_names[0]]

    # Constructs the monolithic function space using the finite element
    # or the set of finite elements, which were previously created

    monolithic_functionSpace = FunctionSpace(mesh_dataClass.mesh, 
    mixed_element)

    # Creates the function for the updated solution, i.e. the vector of 
    # parameters. Then, splits into the individual fields. As the mono-
    # lithic function space has an specific order of fields, retrieves 
    # the individual fields following this order. Does the same for the
    # variations

    trial_functions = TrialFunction(monolithic_functionSpace)

    monolithic_solution = Function(monolithic_functionSpace)

    solution_functions = []

    variation_functions = []

    if len(fields_names)>1:

        solution_functions = split(monolithic_solution)

        variation_functions = split(TestFunction(monolithic_functionSpace
        ))

    else:

        solution_functions = [monolithic_solution]

        variation_functions = [TestFunction(monolithic_functionSpace)]

    # Organizes them into dictionaries

    solution_fields = dict()

    variation_fields = dict()

    fields_namesDict = dict()

    for i in range(len(fields_names)):

        field_name = fields_names[i]

        # Retrives the individual field and its variation

        solution_fields[field_name] = solution_functions[i]

        variation_fields[field_name] = variation_functions[i]

        fields_namesDict[field_name] = i

    if verbose:

        print("Finishes creating the function spaces\n")

    return (monolithic_functionSpace, monolithic_solution, fields_names, 
    solution_fields, variation_fields, trial_functions, fields_namesDict)

# Defines a function to construct a dictionary of elements from a dic-
# tionary of finite elements' instructions

def construct_elementsDictionary(elements_dictionary, mesh_dataClass):

    # Sets a list of allowed interpolation functions

    allowed_interpolationFunction = ["CG", "DG", "Lagrange"]

    # Initializes a list of names of the fields

    fields_names = []

    # Iterates through the dictionary of elements
    
    for field_name, element_dictionary in elements_dictionary.items():

        # Verifies if this is already a fenics element. Tests by veri-
        # fying if has a cell attribute

        try:

            cell = element_dictionary.cell()

            # Updates the list of fields' names

            fields_names.append(field_name)

        # If not, it must be a dictionary with finite element's informa-
        # tion

        except:

            # Verifies if the element is in fact a function space

            if isinstance(element_dictionary, FunctionSpace):

                # Gets the element from the function space

                elements_dictionary[field_name] = (
                element_dictionary.ufl_element())

            # Otherwise, creates the element

            else:

                if not isinstance(element_dictionary, dict):

                    raise TypeError("The dictionary to construct the f"+
                    "inite element of the field '"+str(field_name)+"' "
                    "is not a dictionary. Check out: "+str(
                    element_dictionary))
                
                # Gets the polynomial degree. Uses second degree polyno-
                # mial as default

                polynomial_degree = 2

                if "polynomial degree" in element_dictionary:

                    polynomial_degree = element_dictionary["polynomial"+
                    " degree"]
                
                # Gets the interpolation function

                interpolation_function = ""

                if "interpolation function" in element_dictionary:

                    interpolation_function = element_dictionary["inter"+
                    "polation function"]

                    # Verifies if this interpolation function is one of 
                    # the admissible ones

                    if not (interpolation_function in (
                    allowed_interpolationFunction)):
                        
                        raise NameError("The interpolation function '"+
                        str(interpolation_function)+"' is not one of t"+
                        "he admissible functions. Check the list ahead"+
                        ": "+str(allowed_interpolationFunction))
                    
                elif "shape function" in element_dictionary:

                    interpolation_function = element_dictionary["shape"+
                    " function"]

                    # Verifies if this interpolation function is one of 
                    # the admissible ones

                    if not (interpolation_function in (
                    allowed_interpolationFunction)):
                        
                        raise NameError("The interpolation function '"+
                        str(interpolation_function)+"' is not one of t"+
                        "he admissible functions. Check the list ahead"+
                        ": "+str(allowed_interpolationFunction))
                    
                else:

                    # Uses continuous Galerkin as default

                    interpolation_function = "CG"

                # Gets the element type (scalar, vector, tensor...)

                if "field type" in element_dictionary:

                    element_type = element_dictionary["field type"]

                    # If the field is a scalar-valued field

                    if element_type=="scalar":

                        # Creates the element

                        elements_dictionary[field_name] = FiniteElement(
                        interpolation_function, 
                        mesh_dataClass.mesh.ufl_cell(), polynomial_degree)

                    # If the field is a vector-valued field

                    elif element_type=="vector":

                        # Creates the element

                        elements_dictionary[field_name] = VectorElement(
                        interpolation_function, 
                        mesh_dataClass.mesh.ufl_cell(), polynomial_degree)

                    # If the field is a tensor-valued field

                    elif element_type=="tensor":

                        # Creates the element

                        elements_dictionary[field_name] = TensorElement(
                        interpolation_function, 
                        mesh_dataClass.mesh.ufl_cell(), polynomial_degree)

                    else:

                        raise NameError("There is no '"+str(element_type
                        )+"' field implemented. A field can be either "+
                        "'scalar', 'vector', or 'tensor' for the eleme"+
                        "nt to be dynamically generated")
                    
                else:

                    raise KeyError("There is no key 'field type' in th"+
                    "e dictionary to create finite elements.")
                
            # Updates the list of fields' names

            fields_names.append(field_name)
        
    # Returns the finite elements dictionary

    return elements_dictionary, fields_names

########################################################################
#                        Time stepping classes                         #
########################################################################

# Defines an abstract class as template for the time-stepping algorithms

class TimeSteppingClasses(ABC):

    @abstractmethod

    # The following methods have the pass argument only because they 
    # will be defined in the child classes

    def update(self, time_value):

        pass

# Defines a class to store macroquantities as constants, then, it is up-
# dated during the time-stepping algorithm. This class retrieves the da-
# ta from txt files, convert them to lists and, then, uses them to as-
# sign to FEniCS constants

class MacroQuantitiesInTime(TimeSteppingClasses):

    def __init__(self, dictionary_VariablesToFiles, time_tolerance=1E-5):

        # Stores the dictionary of variables

        self.variables = dictionary_VariablesToFiles

        # Initializes a list of time keys. All the macro quantities must
        # have the same time points

        self.time_keys = []

        # Iterates through the dictionary of file names to read them in-
        # to lists

        for variable, file_name in self.variables.items():

            # Converts the lists into dictionaries

            self.variables[variable] = (file_tools.list_toDict(
            file_tools.txt_toList(file_name)))

            # Updates the time keys

            if len(self.time_keys)==0:

                self.time_keys = list(self.variables[variable].keys())

            # Otherwise, checks if the time keys are equal given a tole-
            # rance

            else:

                test_timeKeysConsistency(self.time_keys, list(
                self.variables[variable].keys()), time_tolerance=
                time_tolerance, variable_name=variable)

            # Initializes the variable in the class data structure to u-
            # se it after as class.variable. Uses FEniCS constant, so it
            # might be just assigned later to spare compilation time

            setattr(self, variable, Constant(self.variables[variable][
            self.time_keys[0]]))

    # Defines a function to update the macro quantitites given the cur-
    # rent time value

    def update(self, time_value):
        
        # Iterates through the variables

        for variable, variable_dataDict in self.variables.items():

            # Gets the keys of the variable data dictionary

            data_keys = list(variable_dataDict.keys())

            # Gets the closest key to the given time value

            time_key = numerical_tools.closest_numberInList(time_value,
            data_keys)

            # Gets the data value for this time step to the variable 
            # contained in this class with this variable's name. Uses 
            # the assign method, since this variables has already been
            # created as a fenics Constant

            getattr(self, variable).assign(Constant(variable_dataDict[
            time_key]))

            # Deletes this key-data pair from the data dictionary

            variable_dataDict.pop(time_key)

# Defines a function to test if the time keys are the same comparing the
# standard with a given set

def test_timeKeysConsistency(standard_timeKeys, compared_timeKeys, 
variable_name="NoGivenName", time_tolerance=1E-5):

    # Checks if it has the same number of points as the first
    # one

    if len(compared_timeKeys)!=len(standard_timeKeys):

        raise ValueError("All the macro quantities must have the same "+
        "number of time points. The first macro quantity had "+str(len(
        standard_timeKeys))+" time points, whereas quantity "+str(
        variable_name)+" has "+str(len(compared_timeKeys)+" time points"))
    
    # If they have the same number of time points, checks if
    # they are the same points

    for i in range(len(standard_timeKeys)):

        if (abs(standard_timeKeys[i]-compared_timeKeys[i])>
        time_tolerance):
            
            raise ValueError("The time points must be equal for all ma"+
            "cro quantities. The "+str(i)+"-th point of the first quan"+
            "tity is "+str(standard_timeKeys[i])+" whereas the same po"+
            "int in the "+str(variable_name)+" macro quantity is "+str(
            compared_timeKeys[i]))

# Defines a function to transform a key to ascii pattern

def convert_stringToASCII(string_toConvert):

    # Initializes the converted string

    converted_string = ""

    # Iterates through the characters

    for character in string_toConvert:

        # Checks if it is a space

        if character==" " or character==' ':

            # Substitutes the blank space for an underline

            converted_string += "_"

        # Checks if it is ascii. If the character encoding is larger
        # than 128, just ignore it, for it's not ascii

        elif ord(character)<128:

            converted_string += character

    # If the converted string is empty, returns a basic string

    if len(converted_string)==0:

        return "foo_variable"

    # Returns the string

    return converted_string
        
########################################################################
#                       Saving of stress measures                      #
########################################################################

# Defines a function to get, project and save a stress field

def save_stressField(output_object, field, time, flag_parentMeshReuse,
stress_solutionPlotNames, stress_name, stress_method, fields_namesDict):

    # If the flag to reuse parent mesh information is true, just save 
    # the information given in the output object

    if flag_parentMeshReuse:

        output_object.parent_toChildMeshResult.rename(
        *stress_solutionPlotNames)

        output_object.result.write(
        output_object.parent_toChildMeshResult, time)

        return output_object
    
    # Verifies if the output object has the attribute with the names of 
    # the required fields

    if not hasattr(output_object, "required_fieldsNames"):

        raise AttributeError("The class of data for the post-process o"+
        "f saving the stress field does not have the attribute 'requir"+
        "ed_fieldsNames'. This class must have it")

    # Verifies if the domain is homogeneous

    if isinstance(output_object.constitutive_model, dict):

        # Initializes a list of pairs of constitutive models and inte-
        # gration domain

        integration_pairs = []

        # If the domain is heterogeneous, the stress field must be pro-
        # jected for each subdomain

        for subdomain, local_constitutiveModel in output_object.constitutive_model.items():

            # Gets the fields for this constitutive model

            retrieved_fields = select_fields(field, 
            output_object.required_fieldsNames[subdomain], 
            fields_namesDict)

            # Gets the stress field

            stress_field = programming_tools.get_result(getattr(
            local_constitutiveModel, stress_method)(retrieved_fields), 
            stress_name)

            # Verifies if more than one physical group is given for the
            # same constitutive model

            if isinstance(subdomain, tuple):

                # Iterates though the elements of the tuple

                for sub in subdomain:

                    # Converts the subdomain to an integer tag

                    sub = variational_tools.verify_physicalGroups(sub, 
                    output_object.physical_groupsList, 
                    output_object.physical_groupsNamesToTags,
                    throw_error=False)

                    # Checks if this subdomain is in the domain physical
                    # groups 

                    if sub in output_object.physical_groupsList:

                        # Adds this pair of constitutive model and inte-
                        # gration domain to the list of such pairs

                        integration_pairs.append([stress_field, sub])

            else:

                # Converts the subdomain to an integer tag

                subdomain = variational_tools.verify_physicalGroups(
                subdomain, output_object.physical_groupsList, 
                output_object.physical_groupsNamesToTags, throw_error=
                False)

                # Checks if this subdomain is in the domain physical
                # groups 

                if subdomain in output_object.physical_groupsList:

                    # Adds this pair of constitutive model and integra-
                    # tion domain to the list of such pairs

                    integration_pairs.append([stress_field, subdomain])

        # Projects this piecewise continuous field of stress into a FE 
        # space

        stress_fieldFunction = variational_tools.project_piecewiseField(
        integration_pairs, output_object.dx, output_object.W, 
        output_object.physical_groupsList, 
        output_object.physical_groupsNamesToTags, solution_names=
        stress_solutionPlotNames)

        # Saves the field into the sharable result with a submesh

        output_object.parent_toChildMeshResult = stress_fieldFunction

        # Writes the field to the file

        output_object.result.write(stress_fieldFunction, time)

    else:

        # Gets the fields for this constitutive model

        retrieved_fields = select_fields(field, 
        output_object.required_fieldsNames, fields_namesDict)

        # Gets the stress field

        stress_field = programming_tools.get_result(getattr(
        output_object.constitutive_model, stress_method)(field), 
        stress_name)

        # Projects the stress into a function

        stress_fieldFunction = project(stress_field, output_object.W)

        stress_fieldFunction.rename(*stress_solutionPlotNames)

        # Saves the field into the sharable result with a submesh

        output_object.parent_toChildMeshResult = stress_fieldFunction

        # Writes the field to the file

        output_object.result.write(stress_fieldFunction, time)

    return output_object

########################################################################
#     Creation and pre-evaluation of discontinuous function spaces     #
########################################################################

# Defines a function to create a R^3 -> R discontinuous function space
# and populates it with values using physical groups information. The 
# keys are the physical groups tags and the values are the values of the
# property to be interpolated using a function space

def physical_groupToDGSpace(physical_dictionary, mesh, mesh_function):

    # Creates a discontinuous Galerkin function space

    DG_functionSpace = FunctionSpace(mesh, "DG", 0)

    # Creates the function with parameters vector

    DG_function = Function(DG_functionSpace)

    # Iterates through the elements of the mesh

    for cell in cells(mesh):

        # Gets the physical group tag (numerical) to which this element 
        # belongs

        physical_groupTag = mesh_function[cell]

        # Iterates through the keys (physical groups numerical tags) and
        # the values of the physical dictionary

        #flag_test = True

        try:

            DG_function.vector()[cell.index()] = physical_dictionary[
            physical_groupTag]
            
        except:

            raise ValueError("The dictionary of properties "+str(
            physical_dictionary)+" does not contain the physical group"+
            " tag "+str(physical_groupTag))
            
    # Returns the function

    return DG_function