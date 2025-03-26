# Routine to store some methods to post-process solution inside the 
# pseudotime stepping methods

from dolfin import *

import copy

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.variational_tools as variational_tools

import source.tool_box.homogenization_tools as homogenization_tools

########################################################################
#                   Post-processing tools selection                    #
########################################################################

# Defines a function to select the tools of post-process. The variable
# post_processes is a dictionary, where the keys are the name of the 
# processes, and the values are dictionaries of additional information.
# The dictionary of additional information has the name of each additio-
# nal information as keys and the informations themselves as values

def post_processingSelectionSingleField(post_processes, mesh,
constitutive_model, dx):

    # Constructs a dictionary of post processess that are implemented. 
    # If you implement any new post-processing tool, you have to append 
    # it to this list

    available_processes = dict()
    
    # For each post-process (called by its key), the initializer must be
    # supplied, then the method itself, and there is a set of additional 
    # information that might be required, such as file names and other
    # questions. The fourth component is a list of variables that are 
    # supplied directly by the code, such as mesh
    
    available_processes["save field"] = [initialize_fieldSaving, 
    update_fieldSaving, ["directory path", "file name"], []] 
    
    available_processes["save stress field"] = [
    initialize_cauchyStressSaving, update_cauchyStressSaving, ["direct"+
    "ory path", "file name", "polynomial degree"], [mesh, 
    constitutive_model, dx]]
    
    available_processes["homogenize field"] = [
    initialize_fieldHomogenization, update_fieldHomogenization, ["dir"+
    "ectory path", "file name", "subdomain"], [dx]]
    
    available_processes["homogenize field's gradient"] = [
    initialize_gradientFieldHomogenization, 
    update_gradientFieldHomogenization, ["directory path", "file name", 
    "subdomain"], [dx]]

    # Gets the names of the available processes
   
    names_availableProcesses = available_processes.keys()

    # Iterates through the dictionary of wanted processes
    
    for process_name, additional_informationDict in post_processes.items():

        # Verifies if this process is implemented
        
        if not (process_name in names_availableProcesses):

            raise NameError("The process '"+str(process_name)+"' is no"+
            " timplemented. The name might be wrong or it might not be"+
            " implemented. Check out the list of available post-proces"+
            "sing methods:\n"+str(list(names_availableProcesses)))

        # Initializes a list of additional information or data to be 
        # passed to the method
        
        method_data = [info for info in available_processes[process_name
        ][2]]
        
        # Iterates through the possible additional information names

        for i in range(len(available_processes[process_name][2])):

            # Initializes a flag to inform if this information exists in
            # the given information

            info_existence = False 

            # Iterates through the additional information dictionary

            for additional_info, info in additional_informationDict.items():

                if available_processes[process_name][2][i]==additional_info:

                    # Updates the existence flag and substitute the 
                    # method data for the actual data

                    info_existence = True 

                    method_data[i] = copy.deepcopy(info)

                    # Deletes this key

                    del additional_informationDict[additional_info]

                    break 

            # If this information has not been found

            if not info_existence:

                raise NameError("The additional information '"+str(
                available_processes[process_name][2][i])+"' has not be"+
                "en found for the "+str(process_name)+" process.")
            
        # Verifies if any keys have been left out in the original dic-
        # tionary of additional information

        if len(list(additional_informationDict.keys()))>0:

            raise KeyError("The additional infos "+str(list(
            additional_informationDict.keys()))+" are not valid additi"+
            "onal information keys for the '"+str(process_name)+"' pro"+
            "cess.")
        
        # Substitute the true data into the additional information list,
        # and adds the initialization and update methods

        post_processes[process_name] = [available_processes[process_name
        ][0], available_processes[process_name][1], copy.deepcopy(
        method_data), available_processes[process_name][3]]

    # Returns the new and complete dictionary of processes

    return post_processes

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

def update_fieldSaving(file, field, time):

    # Writes the field to the file

    file.write(field, time)

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

    output_object = [file, W, constitutive_model, dx]

    return output_object

# Defines a function to update the Cauchy stress field

def update_cauchyStressSaving(output_object, field, time):

    # Recovers the object items

    file, W, constitutive_model, dx = output_object

    # Verifies if the domain is homogeneous

    if isinstance(constitutive_model, dict):

        # Initializes the Cauchy stress function

        cauchy_stressFunction = Function(W)

        # If the domain is heterogeneous, the stress field must be pro-
        # jected for each subdomain

        for subdomain, local_constitutiveModel in constitutive_model.items():

            # Gets the Cauchy stress field

            cauchy_stress = local_constitutiveModel.cauchy_stress(field)

            # Projects the cauchy stress into a function

            cauchy_stressProjected = variational_tools.projection_overRegion(
            cauchy_stress, W, dx, subdomain)

            # Updates the parameters vector of the FEM interpolation of 
            # the Cauchy stress 

            cauchy_stressFunction.vector()[:] = (
            cauchy_stressFunction.vector()[:]+cauchy_stressProjected.vector()[:])

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

def update_fieldHomogenization(output_object, field, time):

    # Gets the data

    homogenized_fieldList = output_object[0]

    inverse_volume = output_object[1]

    dx = output_object[2]
    
    subdomain = output_object[3]

    file_name = output_object[4]

    # Homogenizes the field and updates the list of homogenized field
    # along time

    homogenized_fieldList = homogenization_tools.homogenize_genericField(
    field, homogenized_fieldList, time, inverse_volume, dx, subdomain, 
    file_name)
    
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

def update_gradientFieldHomogenization(output_object, field, time):

    # Gets the gradient of the field

    grad_field = grad(field)

    # Gets the homogenization of the gradient

    output_object = update_fieldHomogenization(output_object, 
    grad_field, time)

    return output_object