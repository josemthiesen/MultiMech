# Routine to store functions to help with functional analysis tasks

from dolfin import *

########################################################################
#     Creation and pre-evaluation of discontinuous function spaces

# Defines a function to create a R^3 -> R discontinuous function space
# and populates it with values using physical groups information

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

        for key, value in physical_dictionary.items():

            # If the key matches the physical group tag, the value is a-
            # located into the function space

            if key==physical_groupTag:

                DG_function.vector()[cell.index()] = value
            
    # Returns the function

    return DG_function