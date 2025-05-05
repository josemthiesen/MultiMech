# Routine to store functions to help with functional analysis tasks and
# programmatic tools

from dolfin import *

from abc import ABC, abstractmethod

import source.tool_box.file_handling_tools as file_tools

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

    def __init__(self, dictionary_VariablesToFiles, ):

        # Stores the dictionary of variables

        self.variables = dictionary_VariablesToFiles

        # Iterates through the dictionary of file names to read them in-
        # to lists

        for variable, file_name in self.variables.items():

            # Converts the lists into dictionaries

            self.variables[variable] = (file_tools.list_toDict(
            file_tools.txt_toList(file_name)))

            # Initializes the variable in the class data structure to u-
            # se it after as class.variable

            setattr(self, variable, self.variables[variable][list(
            self.variables[variable].keys())[0]])

    # Defines a function to update the macro quantitites given the cur-
    # rent time value

    def update(self, time_value):
        
        # Iterates through the variables

        for variable, variable_dataDict in self.variables.items():

            pass

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