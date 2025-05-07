# Routine to store functions to help with functional analysis tasks and
# programmatic tools

from dolfin import *

from abc import ABC, abstractmethod

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.numerical_tools as numerical_tools

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