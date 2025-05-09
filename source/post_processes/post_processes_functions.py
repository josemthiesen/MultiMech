# Routine to store some methods to post-process solution inside the 
# pseudotime stepping methods

from dolfin import *

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.homogenization_tools as homogenization_tools

import source.tool_box.functional_tools as functional_tools

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
    
    # Assembles the output. This post-process does not have a variable
    # that can be shared with a submesh

    class OutputObject:

        def __init__(self, file):
            
            self.result = file

    # Initializes the file

    file = XDMFFile(file_name)

    output_object = OutputObject(file)

    return output_object

# Defines a function to update the file with the field

def update_fieldSaving(output_object, field, field_number, time):
    
    print("Updates the saving of the "+str(field_number)+" field\n")

    # If the problem has a single field

    if field_number==-1:

        # Writes the field to the file

        output_object.result.write(field, time)

        return output_object
    
    # If the problem has multiple fields

    else:

        # Writes the field to the file

        output_object.result.write(field[field_number], time)

        return output_object

########################################################################
#                          Cauchy stress field                         #
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

    # Assembles the file and the function space into a class. This post-
    # process does have a variable that can be shared with a submesh, 
    # and it is the stress field

    class OutputObject:

        def __init__(self, file, W, constitutive_model, dx, 
        physical_groupsList, physical_groupsNamesToTags, 
        parent_toChildMeshResult):
            
            self.result = file 

            self.W = W 

            self.constitutive_model = constitutive_model

            self.dx = dx 

            self.physical_groupsList = physical_groupsList 

            self.physical_groupsNamesToTags = physical_groupsNamesToTags
            
            # Defines a sharable result between a parent mesh and a sub-
            # mesh

            self.parent_toChildMeshResult = parent_toChildMeshResult

    output_object = OutputObject(file, W, constitutive_model, dx, 
    physical_groupsList, physical_groupsNamesToTags, 0.0)

    return output_object

# Defines a function to update the Cauchy stress field

def update_cauchyStressSaving(output_object, field, field_number, time, 
flag_parentMeshReuse=False):
    
    print("Updates the saving of the Cauchy stress field\n")
    
    return functional_tools.save_stressField(output_object, field, time, 
    flag_parentMeshReuse, ["Cauchy stress", "stress"], "cauchy", "cauc"+
    "hy_stress")

########################################################################
#                      Couple Cauchy stress field                      #
########################################################################

# Defines a function to initialize the couple Cauchy stress field file

def initialize_coupleCauchyStressSaving(data, direct_codeData, 
submesh_flag):

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

    # Assembles the file and the function space into a class. This post-
    # process does have a variable that can be shared with a submesh, 
    # and it is the couple stress field

    class OutputObject:

        def __init__(self, file, W, constitutive_model, dx, 
        physical_groupsList, physical_groupsNamesToTags,
        parent_toChildMeshResult):
            
            self.result = file 

            self.W = W 

            self.constitutive_model = constitutive_model

            self.dx = dx 

            self.physical_groupsList = physical_groupsList 

            self.physical_groupsNamesToTags = physical_groupsNamesToTags
            
            # Defines a sharable result between a parent mesh and a sub-
            # mesh

            self.parent_toChildMeshResult = parent_toChildMeshResult

    output_object = OutputObject(file, W, constitutive_model, dx, 
    physical_groupsList, physical_groupsNamesToTags, 0.0)

    return output_object

# Defines a function to update the couple Cauchy stress field

def update_coupleCauchyStressSaving(output_object, field, field_number, 
time, flag_parentMeshReuse=False):
    
    print("Updates the saving of the couple Cauchy stress field\n")

    return functional_tools.save_stressField(output_object, field, time, 
    flag_parentMeshReuse, ["Couple Cauchy stress", "stress"], "couple_"+
    "cauchy", "cauchy_stress")

########################################################################
#                  First Piola-Kirchhoff stress field                  #
########################################################################

# Defines a function to initialize the first Piola-Kirchhoff stress 
# field file

def initialize_firstPiolaStressSaving(data, direct_codeData, 
submesh_flag):

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

    # Assembles the file and the function space into a class. This post-
    # process does have a variable that can be shared with a submesh, 
    # and it is the stress field

    class OutputObject:

        def __init__(self, file, W, constitutive_model, dx, 
        physical_groupsList, physical_groupsNamesToTags, 
        parent_toChildMeshResult):
            
            self.result = file 

            self.W = W 

            self.constitutive_model = constitutive_model

            self.dx = dx 

            self.physical_groupsList = physical_groupsList 

            self.physical_groupsNamesToTags = physical_groupsNamesToTags
            
            # Defines a sharable result between a parent mesh and a sub-
            # mesh

            self.parent_toChildMeshResult = parent_toChildMeshResult

    output_object = OutputObject(file, W, constitutive_model, dx, 
    physical_groupsList, physical_groupsNamesToTags, 0.0)

    return output_object

# Defines a function to update the first Piola-Kirchhoff stress field

def update_firstPiolaStressSaving(output_object, field, field_number, 
time, flag_parentMeshReuse=False):
    
    print("Updates the saving of the first Piola-Kirchhoff stress fiel"+
    "d\n")
    
    return functional_tools.save_stressField(output_object, field, time, 
    flag_parentMeshReuse, ["First Piola-Kirchhoff stress", "stress"], 
    "first_piola_kirchhoff", "first_piolaStress")

########################################################################
#              Couple first Piola-Kirchhoff stress field               #
########################################################################

# Defines a function to initialize the couple first Piola-Kirchhoff 
# stress field file

def initialize_coupleFirstPiolaStressSaving(data, direct_codeData, 
submesh_flag):

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

    # Assembles the file and the function space into a class. This post-
    # process does have a variable that can be shared with a submesh, 
    # and it is the stress field

    class OutputObject:

        def __init__(self, file, W, constitutive_model, dx, 
        physical_groupsList, physical_groupsNamesToTags, 
        parent_toChildMeshResult):
            
            self.result = file 

            self.W = W 

            self.constitutive_model = constitutive_model

            self.dx = dx 

            self.physical_groupsList = physical_groupsList 

            self.physical_groupsNamesToTags = physical_groupsNamesToTags
            
            # Defines a sharable result between a parent mesh and a sub-
            # mesh

            self.parent_toChildMeshResult = parent_toChildMeshResult

    output_object = OutputObject(file, W, constitutive_model, dx, 
    physical_groupsList, physical_groupsNamesToTags, 0.0)

    return output_object

# Defines a function to update the couple first Piola-Kirchhoff stress 
# field

def update_coupleFirstPiolaStressSaving(output_object, field, 
field_number, time, flag_parentMeshReuse=False):
    
    print("Updates the saving of the couple first Piola-Kirchhoff stre"+
    "ss field\n")
    
    return functional_tools.save_stressField(output_object, field, time, 
    flag_parentMeshReuse, ["Couple first Piola-Kirchhoff stress", "str"+
    "ess"], "couple_first_piola_kirchhoff", "first_piolaStress")

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

    physical_groupsList = direct_codeData[1] 
    
    physical_groupsNamesToTags = direct_codeData[2]

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

        # Checks for any elements being strings

        for i in range(len(subdomain)):
            
            if isinstance(subdomain[i], str):

                subdomain[i] = variational_tools.verify_physicalGroups(
                subdomain[i], physical_groupsList, 
                physical_groupsNamesToTags=physical_groupsNamesToTags)

        # Converts to tuple

        subdomain = tuple(subdomain)

    if isinstance(subdomain, int):

        volume = assemble(1*dx(subdomain))

    elif isinstance(subdomain, tuple):

        for sub in subdomain:

            volume += assemble(1*dx(sub))

    # Otherwise, integrates over the whole domain to get the volume

    elif isinstance(subdomain, str):

        if len(subdomain)==0:

            volume = assemble(1*dx)

        else:

            subdomain = variational_tools.verify_physicalGroups(
            subdomain, physical_groupsList, physical_groupsNamesToTags=
            physical_groupsNamesToTags)

            volume = assemble(1*dx(subdomain))

    # Initializes the homogenized field list

    homogenized_fieldList = []

    # Gets the name of the file with the path to it

    file_name = file_tools.verify_path(parent_path, file_name)

    # Assembles the output. This post-process does not have a variable
    # that can be shared with a submesh


    class OutputObject:

        def __init__(self, homogenized_fieldList, inverse_volume, dx, 
        subdomain, file_name):
            
            self.result = homogenized_fieldList

            self.inverse_volume = inverse_volume

            self.dx = dx 

            self.subdomain = subdomain 

            self.file_name = file_name

    output_object = OutputObject(homogenized_fieldList, (1.0/volume), dx, 
    subdomain, file_name)

    return output_object

# Defines a function to update the homogenized field

def update_fieldHomogenization(output_object, field, field_number, time):
    
    print("Updates the homogenization of the "+str(field_number)+" fie"+
    "ld\n")

    # If the problem has a single field

    if field_number==-1:

        # Homogenizes the field and updates the list of homogenized field
        # along time

        output_object.result = homogenization_tools.homogenize_genericField(
        field, output_object.result, time, output_object.inverse_volume, 
        output_object.dx, output_object.subdomain, 
        output_object.file_name)

        return output_object

    # If the problem has multiple fields

    else:

        # Homogenizes the field and updates the list of homogenized field
        # along time

        output_object.result = homogenization_tools.homogenize_genericField(
        field[field_number], output_object.result, time, 
        output_object.inverse_volume, output_object.dx, 
        output_object.subdomain, output_object.file_name)

        return output_object

# Defines a function to initialize the homogenization of the gradient of 
# a field

def initialize_gradientFieldHomogenization(data, direct_codeData, 
submesh_flag):

    output = initialize_fieldHomogenization(data, direct_codeData, 
    submesh_flag)

    return output

# Defines a function to update the homogenization of the gradient of a 
# field

def update_gradientFieldHomogenization(output_object, field, 
field_number, time):
    
    print("Updates the homogenization of the gradient of the "+str(
    field_number)+" field\n")

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

# Defines a function to initialize the homogenized value of the first 
# Piola-Kirchhof

def initialize_firstPiolaHomogenization(data, direct_codeData, 
submesh_flag):

    # Gets the directory and the name of the file

    parent_path = data[0]

    file_name = data[1]

    # Gets the subdomain to integrate

    subdomain = data[2]

    # Gets the integration measure

    dx = direct_codeData[0]

    physical_groupsList = direct_codeData[1] 
    
    physical_groupsNamesToTags = direct_codeData[2]

    constitutive_model = direct_codeData[3]

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

    if isinstance(subdomain, int):

        volume = assemble(1*dx(subdomain))

    elif isinstance(subdomain, tuple):

        for sub in subdomain:

            if isinstance(sub, str):

                volume += assemble(1*dx(variational_tools.verify_physicalGroups(
                sub, physical_groupsList, physical_groupsNamesToTags=
                physical_groupsNamesToTags)))

            else:

                volume += assemble(1*dx(sub))

    # Otherwise, integrates over the whole domain to get the volume

    elif isinstance(subdomain, str):

        if len(subdomain)==0:

            volume = assemble(1*dx)

        else:

            volume = assemble(1*dx(variational_tools.verify_physicalGroups(
            subdomain, physical_groupsList, physical_groupsNamesToTags=
            physical_groupsNamesToTags)))

    # Initializes the homogenized field list

    homogenized_firstPiolaList = []

    # Gets the name of the file with the path to it

    file_name = file_tools.verify_path(parent_path, file_name)

    # Assembles the output. This post-process does not have a variable
    # that can be shared with a submesh

    class OutputObject:

        def __init__(self, homogenized_firstPiolaList, inverse_volume, dx, 
        subdomain, file_name, constitutive_model, physical_groupsList,
        physical_groupsNamesToTags):
            
            self.result = homogenized_firstPiolaList

            self.inverse_volume = inverse_volume

            self.dx = dx 

            self.subdomain = subdomain 

            self.file_name = file_name

            self.constitutive_model = constitutive_model

            self.physical_groupsList = physical_groupsList
    
            self.physical_groupsNamesToTags = physical_groupsNamesToTags

    output_object = OutputObject(homogenized_firstPiolaList, (1.0/volume
    ), dx, subdomain, file_name, constitutive_model, physical_groupsList,
    physical_groupsNamesToTags)

    return output_object

# Defines a function to update the homogenization of the first Piola-
# Kirchhof

def update_firstPiolaHomogenization(output_object, field, field_number, 
time):
    
    print("Updates the homogenization of the first Piola-Kirchhoff str"+
    "ess field\n")

    output_object.result = homogenization_tools.homogenize_stressTensor(
    field, output_object.constitutive_model, "first_piola_kirchhoff", 
    "first_piolaStress", output_object.result, time, 
    output_object.inverse_volume, output_object.dx, 
    output_object.subdomain,output_object.file_name, 
    output_object.physical_groupsList, 
    output_object.physical_groupsNamesToTags)

    return output_object

# Defines a function to initialize the homogenized value of the first 
# couple Piola-Kirchhof

def initialize_coupleFirstPiolaHomogenization(data, direct_codeData, 
submesh_flag):
    
    return initialize_firstPiolaHomogenization(data, direct_codeData, 
    submesh_flag)

# Defines a function to update the homogenization of the couple first 
# Piola-Kirchhof

def update_coupleFirstPiolaHomogenization(output_object, field, 
field_number, time):
    
    print("Updates the homogenization of the couple first Piola-Kirchh"+
    "off stress field\n")

    output_object.result = homogenization_tools.homogenize_stressTensor(
    field, output_object.constitutive_model, "couple_first_piola_kirch"+
    "hoff", "first_piolaStress", output_object.result, time, 
    output_object.inverse_volume, output_object.dx, 
    output_object.subdomain,output_object.file_name, 
    output_object.physical_groupsList, 
    output_object.physical_groupsNamesToTags)

    return output_object