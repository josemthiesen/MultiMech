# Routine to store some methods to post-process solution inside the 
# pseudotime stepping methods

from dolfin import *

import copy

import source.tool_box.file_handling_tools as file_tools

########################################################################
#                   Post-processing tools selection                    #
########################################################################

# Defines a function to select the tools of post-process. The variable
# post_processes is a dictionary, where the keys are the name of the 
# processes, and the values are dictionaries of additional information.
# The dictionary of additional information has the name of each additio-
# nal information as keys and the informations themselves as values

def post_processingSelectionSingleField(post_processes):

    # Constructs a dictionary of post processess that are implemented. 
    # If you implement any new post-processing tool, you have to append 
    # it to this list

    available_processes = dict()
    
    # For each post-process (called by its key), the initializer must be
    # supplied, then the method itself, and there is a set of additional 
    # information that might be required, such as file names and other
    # questions
    
    available_processes["save field"] = [initialize_fieldSaving, 
    update_fieldSaving, ["directory path", "file name"]] 
    
    """available_processes["save dual field"] = [["directory path", "file"+
    " name", "constitutive tensor"]]
    
    available_processes["homogenize field"] = [["directory path", "fil"+
    "e name"]]
    
    available_processes["homogenize field's gradient"] = [["directory "+
    "path", "file name"]]"""

    # Gets the names of the available processes
   
    names_availableProcesses = available_processes.keys()

    # Iterates through the dictionary of wanted processes
    
    for process_name, additional_informationDict in post_processes.items():

        # Verifies if this process is implemented
        
        if not (process_name in names_availableProcesses):

            raise NameError("The process "+str(process_name)+" is not "+
            "implemented. The name might be wrong or it might not be i"+
            "mplemented. Check out the list of available post-processi"+
            "ng methods:\n"+str(names_availableProcesses))

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

                raise NameError("The additional information "+str(
                available_processes[process_name][2][i])+" has not bee"+
                "n found for the "+str(process_name)+" process.")
            
        # Verifies if any keys have been left out in the original dic-
        # tionary of additional information

        if len(list(additional_informationDict.keys()))>0:

            raise KeyError("The additional infos "+str(list(
            additional_informationDict.keys()))+" are not valid additi"+
            "onal information keys for the "+str(process_name)+" proce"+
            "ss.")
        
        # Substitute the true data into the additional information list,
        # and adds the initialization and update methods

        post_processes[process_name] = [available_processes[process_name
        ][0], available_processes[process_name][1], copy.deepcopy(
        method_data)]

    # Returns the new and complete dictionary of processes

    return post_processes

########################################################################
#                             Field saving                             #
########################################################################

# Defines a function to initialize and save a field

def initialize_fieldSaving(data):

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