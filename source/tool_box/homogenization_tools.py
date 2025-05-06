# Routine to store methods for homogenization

from dolfin import *

import numpy as np

import source.tool_box.file_handling_tools as file_tools

# Defines a function to homogenize a generic field

def homogenize_genericField(field, homogenized_fieldList, time, 
inverse_volume, dx, subdomain, file_name):
    
    # Gets the dimensionality of the field

    dimensionality = 0

    try:

        dimensionality = field.value_shape()

    # If the value shape could not be retrieved, it might be a ufl
    # object, such as a gradient

    except:

        dimensionality = field.ufl_shape

    # Initializes the homogenized value

    homogenized_value = 0.0

    # If the field is a scalar

    if len(dimensionality)==0:

        if isinstance(subdomain, int) or isinstance(subdomain, tuple):

            homogenized_value += inverse_volume*assemble(field*dx(
            subdomain))

        else:

            homogenized_value += inverse_volume*assemble(field*dx)

    # If the field is not scalar, iterates through the dimensions

    else:

        # Makes the homogenized value a list with the format of the 
        # field to be homogenized, e.g. vectors, tensors and so forth

        homogenized_value = np.zeros(dimensionality)

        # Get the possible combinations of indexes

        indexes_combinations = file_tools.get_indexesCombinations(
        dimensionality)

        if isinstance(subdomain, int):

            # Iterates through the combinations

            for index_combination in indexes_combinations:

                homogenized_value[*index_combination] = (inverse_volume*
                assemble(field[*index_combination]*dx(subdomain)))
                
        elif isinstance(subdomain, tuple):

            # Iterates through the combinations

            for index_combination in indexes_combinations:

                for subsubdomain in subdomain:

                    homogenized_value[*index_combination] += (
                    inverse_volume*assemble(field[*index_combination]*
                    dx(subsubdomain)))

        else:

            # Iterates through the combinations

            for index_combination in indexes_combinations:

                homogenized_value[*index_combination] = (inverse_volume*
                assemble(field[*index_combination]*dx))

        # Gets the homogenized value back to a list

        homogenized_value = homogenized_value.tolist()

    # Updates the homogenized field

    homogenized_fieldList.append([time, homogenized_value])

    # Saves the homogenized field to a txt file

    file_tools.list_toTxt(homogenized_fieldList, file_name, 
    add_extension=False)

    return homogenized_fieldList

########################################################################
########################################################################

def get_micro_cell_volume(dx):
    RVE_volume = assemble(1*dx(3)) + assemble(1*dx(4))
    return RVE_volume

def get_homogenized(variable, dx):
    # Calculate the volume of the Representative Volume Element (RVE)
    RVE_volume = get_micro_cell_volume(dx)

    # Initialize the homogenized deformation gradient tensor
    Homogenized_variable = np.zeros(3) 

    for m in range(3): 
        # Compute homogenized components
        Homogenized_variable[m] = (
            (1 / RVE_volume) * assemble(variable[m] * dx(3)) + 
            (1 / RVE_volume) * assemble(variable[m] * dx(4))
        )

    return Homogenized_variable 

def get_homogenized_gradient(grad_variable, dx):
    # Calculate the volume of the Representative Volume Element (RVE)
    RVE_volume = get_micro_cell_volume(dx)

    # Initialize the homogenized deformation gradient tensor
    Homogenized_gradient = np.zeros((3, 3)) 

    # Loop through tensor components
    for m in range(3): 
        for n in range(3):
            # Compute homogenized components
            Homogenized_gradient[m, n] = (
                (1 / RVE_volume) * assemble(grad_variable[m, n] * dx(3)) + 
                (1 / RVE_volume) * assemble(grad_variable[m, n] * dx(4))
            )

    return Homogenized_gradient 

