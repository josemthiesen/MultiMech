# Routine to store functions to help with functional analysis tasks and
# programmatic tools

from dolfin import *

from abc import ABC, abstractmethod

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.numerical_tools as numerical_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.programming_tools as programming_tools

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
        
########################################################################
#                       Saving of stress measures                      #
########################################################################

# Defines a function to get, project and save a stress field

def save_stressField(output_object, field, time, flag_parentMeshReuse,
stress_solutionPlotNames, stress_name, stress_method):

    # If the flag to reuse parent mesh information is true, just save 
    # the information given in the output object

    if flag_parentMeshReuse:

        output_object.parent_toChildMeshResult.rename(
        *stress_solutionPlotNames)

        output_object.result.write(
        output_object.parent_toChildMeshResult, time)

        return output_object

    # Verifies if the domain is homogeneous

    if isinstance(output_object.constitutive_model, dict):

        # Initializes a list of pairs of constitutive models and inte-
        # gration domain

        integration_pairs = []

        # If the domain is heterogeneous, the stress field must be pro-
        # jected for each subdomain

        for subdomain, local_constitutiveModel in output_object.constitutive_model.items():

            # Gets the stress field

            stress_field = programming_tools.get_result(getattr(
            local_constitutiveModel, stress_method)(field), stress_name)

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