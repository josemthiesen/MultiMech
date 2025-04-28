# Routine to store functions to help with functional analysis tasks and
# programmatic tools

from dolfin import *

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