# Routine to store some methods to post-process solution inside the 
# pseudotime stepping methods

from dolfin import *

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.homogenization_tools as homogenization_tools

########################################################################
#                      Post-processing tools list                      #
########################################################################

########################################################################
#                             Field saving                             #
########################################################################

# Defines a function to initialize and save a field

def initialize_fieldSaving(data, direct_codeData, submesh_flag):

    # Gets the directory and the name of the file

    parent_path = data[0]

    file_name = data[1]

    # Gets the name of the file with the path to it

    file_name = file_tools.verify_path(parent_path, file_name)

    # Initializes the file

    file = XDMFFile(file_name)

    return file 

# Defines a function to update the file with the field

def update_fieldSaving(file, field, field_number, time):

    # If the problem has a single field

    if field_number==-1:

        # Writes the field to the file

        file.write(field, time)

        return file
    
    # If the problem has multiple fields

    else:

        # Writes the field to the file

        file.write(field[field_number], time)

        return file

########################################################################
#                             Stress field                             #
########################################################################

# Defines a function to initialize the Cauchy stress field file

def initialize_cauchyStressSaving(data, direct_codeData, submesh_flag):

    # Gets the directory and the name of the file

    parent_path = data[0]

    file_name = data[1]

    # Gets the polynomial degree of the interpolation function

    polynomial_degree = data[2]

    # Gets the mesh, the constitutive model, and the volume integrator
    # from the data directly provided by the code

    mesh = direct_codeData[0]

    constitutive_model = direct_codeData[1]

    dx = direct_codeData[2]

    physical_groupsList = direct_codeData[3] 
    
    physical_groupsNamesToTags = direct_codeData[4]

    # Creates the function space for the stress as a tensor

    W = 0.0

    if polynomial_degree==0:

        W = TensorFunctionSpace(mesh, "DG", 0)

    else:

        W = TensorFunctionSpace(mesh, "CG", polynomial_degree)

    # Gets the name of the file with the path to it

    file_name = file_tools.verify_path(parent_path, file_name)

    # Initializes the file

    file = XDMFFile(file_name)

    # Assembles the file and the function space into a list

    output_object = [file, W, constitutive_model, dx, 
    physical_groupsList, physical_groupsNamesToTags]

    return output_object

# Defines a function to update the Cauchy stress field

def update_cauchyStressSaving(output_object, field, field_number, time):

    # Recovers the object items

    (file, W, constitutive_model, dx, physical_groupsList, 
    physical_groupsNamesToTags) = output_object

    # Verifies if the domain is homogeneous

    if isinstance(constitutive_model, dict):

        # Initializes the Cauchy stress function

        cauchy_stressFunction = Function(W)

        # If the domain is heterogeneous, the stress field must be pro-
        # jected for each subdomain

        for subdomain, local_constitutiveModel in constitutive_model.items():

            # Gets the Cauchy stress field

            cauchy_stress = local_constitutiveModel.cauchy_stress(field)

            # Verifies if more than one physical group is given for the
            # same constitutive model

            if isinstance(subdomain, tuple):

                # Iterates though the elements of the tuple

                for sub in subdomain:

                    # Projects the cauchy stress into a function

                    cauchy_stressProjected = variational_tools.projection_overRegion(
                    cauchy_stress, W, dx, sub, physical_groupsList,
                    physical_groupsNamesToTags)

                    # Updates the parameters vector of the FEM interpolation of 
                    # the Cauchy stress 

                    cauchy_stressFunction.vector()[:] = (
                    cauchy_stressFunction.vector()[:]+
                    cauchy_stressProjected.vector()[:])

            else:

                # Projects the cauchy stress into a function

                cauchy_stressProjected = variational_tools.projection_overRegion(
                cauchy_stress, W, dx, subdomain, physical_groupsList,
                physical_groupsNamesToTags)

                # Updates the parameters vector of the FEM interpolation of 
                # the Cauchy stress 

                cauchy_stressFunction.vector()[:] = (
                cauchy_stressFunction.vector()[:]+
                cauchy_stressProjected.vector()[:])

        # Writes the field to the file

        file.write(cauchy_stressFunction, time)

    else:

        # Gets the Cauchy stress field

        cauchy_stress = constitutive_model.cauchy_stress(field)

        # Projects the cauchy stress into a function

        cauchy_stressFunction = project(cauchy_stress, W)

        # Writes the field to the file

        file.write(cauchy_stressFunction, time)

    return [file, W, constitutive_model, dx]

########################################################################
#                            Homogenization                            #
########################################################################

# Defines a function to initialize the homogenization of the field

def initialize_fieldHomogenization(data, direct_codeData, submesh_flag):

    # Gets the directory and the name of the file

    parent_path = data[0]

    file_name = data[1]

    # Gets the subdomain to integrate

    subdomain = data[2]

    # Gets the integration measure

    dx = direct_codeData[0]

    # Evaluates the volume of the domain

    volume = 0.0

    # If the solution comes from a submesh, there can be no domain

    if submesh_flag:

        if (isinstance(subdomain, int) or isinstance(subdomain, tuple)
        or isinstance(subdomain, list)):
            
            raise ValueError("This solution comes from a submesh and t"+
            "he subdomain "+str(subdomain)+" is solicited. Subdomains "+
            "cannot be used in fields from submeshes, for theses meshe"+
            "s do not have physical groups.")

    # If a physical group of the mesh is given or a tuple of physical 
    # groups

    if isinstance(subdomain, list):

        subdomain = tuple(subdomain)

    if isinstance(subdomain, int) or isinstance(subdomain, tuple):

        volume = assemble(1*dx(subdomain))

    # Otherwise, integrates over the whole domain to get the volume

    else:

        volume = assemble(1*dx)

    # Initializes the homogenized field list

    homogenized_fieldList = []

    # Gets the name of the file with the path to it

    file_name = file_tools.verify_path(parent_path, file_name)

    # Assembles the output

    output = [homogenized_fieldList, (1.0/volume), dx, subdomain, 
    file_name]

    return output

# Defines a function to update the homogenized field

def update_fieldHomogenization(output_object, field, field_number, time):

    # Gets the data

    homogenized_fieldList = output_object[0]

    inverse_volume = output_object[1]

    dx = output_object[2]
    
    subdomain = output_object[3]

    file_name = output_object[4]

    # If the problem has a single field

    if field_number==-1:

        # Homogenizes the field and updates the list of homogenized field
        # along time

        homogenized_fieldList = homogenization_tools.homogenize_genericField(
        field, homogenized_fieldList, time, inverse_volume, dx, 
        subdomain, file_name)
        
        # Assembles the output

        output = [homogenized_fieldList, inverse_volume, dx, subdomain, 
        file_name]

        return output

    # If the problem has multiple fields

    else:

        # Homogenizes the field and updates the list of homogenized field
        # along time

        homogenized_fieldList = homogenization_tools.homogenize_genericField(
        field[field_number], homogenized_fieldList, time, inverse_volume, 
        dx, subdomain, file_name)
        
        # Assembles the output

        output = [homogenized_fieldList, inverse_volume, dx, subdomain, 
        file_name]

        return output

# Defines a function to initialize the homogenization of the gradient of 
# a field

def initialize_gradientFieldHomogenization(data, direct_codeData, 
submesh_flag):

    output = initialize_fieldHomogenization(data, direct_codeData, 
    submesh_flag)

    return output

# Defines a function to update the homogenization of the gradient of a 
# field

def update_gradientFieldHomogenization(output_object, field, field_number, time):

    # Gets the gradient of the field

    grad_field = 0.0

    # If the problem has only one field

    if field_number==-1:

        grad_field = grad(field)

    # If the problem has multiple fields

    else:

        grad_field = grad(field[field_number])

    # Gets the homogenization of the gradient. Sets the field number to
    # -1, for a single field is sent downstream

    output_object = update_fieldHomogenization(output_object, 
    grad_field, -1, time)

    return output_object