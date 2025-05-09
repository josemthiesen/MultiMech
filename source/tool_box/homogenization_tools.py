# Routine to store methods for homogenization

from dolfin import *

import numpy as np

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.programming_tools as programming_tools

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
#                        Stress homogenization                         #
########################################################################

# Defines a function to homogenize a stress tensor given by a constitu-
# tive relation

def homogenize_stressTensor(field, constitutive_model, stress_name,
stress_method, homogenized_firstPiolaList, time, inverse_volume, dx, 
subdomain, file_name, physical_groupsList, physical_groupsNamesToTags):
    
    # Initializes the homogenized tensor

    homogenized_tensor = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 
    0.0]]

    tensor_indexes = [[0,0], [0,1], [0,2], [1,0], [1,1], [1,2], [2,0], [
    2,1], [2,2]]
    
    # Verifies if the domain is homogeneous. If the constitutive model 
    # is a dictionary, the domain in heterogeneous

    if isinstance(constitutive_model, dict):

        # If the domain is heterogeneous, the stress field must be pro-
        # jected for each subdomain

        for subdomain, local_constitutiveModel in constitutive_model.items():

            # Gets the stress tensor field

            stress_field = programming_tools.get_result(
            getattr(local_constitutiveModel, stress_method)(field), 
            stress_name)

            # Verifies if more than one physical group is given for the
            # same constitutive model

            if isinstance(subdomain, tuple):

                # Iterates though the elements of the tuple

                for sub in subdomain:

                    # Checks if the domain is a string

                    if isinstance(sub, str):

                        sub = variational_tools.verify_physicalGroups(
                        sub, physical_groupsList, 
                        physical_groupsNamesToTags=
                        physical_groupsNamesToTags)

                    # Adds the contribution of this domain

                    for index in tensor_indexes:

                        homogenized_tensor[index[0]][index[1]] += (
                        inverse_volume*assemble(stress_field[*index
                        ]*dx(sub)))

            else:

                # Checks if the domain is a string

                if isinstance(subdomain, str):

                    subdomain = variational_tools.verify_physicalGroups(
                    subdomain, physical_groupsList, 
                    physical_groupsNamesToTags= 
                    physical_groupsNamesToTags)

                # Adds the contribution of this domain

                for index in tensor_indexes:

                    homogenized_tensor[index[0]][index[1]] += (
                    inverse_volume*assemble(stress_field[*index]*
                    dx(subdomain)))

    else:

        # Gets the stress tensor field

        stress_field = programming_tools.get_result(
        getattr(local_constitutiveModel, stress_method)(field), 
        stress_name)

        # Adds the contribution of this domain

        for index in tensor_indexes:

            homogenized_tensor[index[0]][index[1]] += (inverse_volume*
            assemble(stress_field[*index]*dx))

    # Adds the homogenized tensor to the list

    homogenized_firstPiolaList.append([time, homogenized_tensor])

    # Saves the homogenized field to a txt file

    file_tools.list_toTxt(homogenized_firstPiolaList, file_name, 
    add_extension=False)

    return homogenized_firstPiolaList