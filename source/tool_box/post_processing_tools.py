# Routine to store some methods to post-process solution inside the 
# pseudotime stepping methods

import copy

import source.post_processes.post_processes_classes as post_classes

import source.tool_box.programming_tools as programming_tools

########################################################################
#                   Post-processing tools selection                    #
########################################################################

# Defines a function to select the tools of post-process. The variable
# post_processes is a dictionary, where the keys are the name of the 
# processes, and the values are dictionaries of additional information.
# The dictionary of additional information has the name of each additio-
# nal information as keys and the informations themselves as values

def post_processingSelectionSingleField(post_processes, context_class):

    # Constructs a dictionary of post processess that are implemented. 
    # If you implement any new post-processing tool, the code will auto-
    # matically retrieve it using the inspect functionality

    available_processes = programming_tools.dispatch_processes(
    post_processes, post_classes, post_classes.PostProcessMethod, 
    reserved_classes=[post_classes.PostProcessMethod, 
    post_classes.PostProcessContext], class_input=context_class)

    # Iterates through the dictionary of wanted processes
    
    for process_name, additional_informationDict in post_processes.items():
        
        # Iterates through the possible additional information names

        for i in range(len(available_processes[process_name
        ].additional_information)):

            # Initializes a flag to inform if this information exists in
            # the given information

            info_existence = False 

            # Iterates through the additional information dictionary

            for additional_info, info in additional_informationDict.items():

                if available_processes[process_name
                ].additional_information[i]==additional_info:

                    # Updates the existence flag and substitute the 
                    # method data for the actual data

                    info_existence = True 

                    available_processes[process_name
                    ].additional_information[i] = copy.deepcopy(
                    info)

                    # Deletes this key

                    del additional_informationDict[additional_info]

                    break 

            # If this information has not been found

            if not info_existence:

                raise NameError("The additional information '"+str(
                available_processes[process_name
                ].additional_information[i])+"' has not been found for"+
                "the "+str(process_name)+" process.")
            
        # Verifies if any keys have been left out in the original dic-
        # tionary of additional information

        if len(list(additional_informationDict.keys()))>0:

            raise KeyError("The additional infos "+str(list(
            additional_informationDict.keys()))+" are not valid additi"+
            "onal information keys for the '"+str(process_name)+"' pro"+
            "cess.")

    # Returns the new and complete dictionary of processes

    return available_processes

# Defines a function to select the tools of post-process for a multi-
# field problem. The variable post_processesUnifield is a list of dicti-
# onaries to work with processes that receive only one field as input, 
# where each component correspond to a field; the keys are the names of 
# the processes, and the values are dictionaries of additional informa-
# tion. A dictionary of fields names as keys and their respective index
# in the solution must be provided, too

def post_processingSelectionMultipleFields(post_processesUnifield,
context_class, fields_names):
    
    # Initializes a list of lists, each sublist has a pair of field num-
    # ber and process name

    post_processesNamesList = []
    
    # Iterates through the number of fields

    for i in range(len(post_processesUnifield)):

        # Gets the number of the field to which this post process be-
        # longs

        try:

            field_number = fields_names[copy.deepcopy(
            post_processesUnifield[i][0])]

        except IndexError:

            raise IndexError("The post-processes list for multiple fie"+
            "ld physics must be a list of lists, i.e. each process is "+
            "a sublist with the first component being the number or th"+
            "e name of the field that is used for this particular post"+
            "-process and the second value is the dictionary of the pr"+
            "ocess")
        
        except KeyError:

            if isinstance(post_processesUnifield[i][0], int):

                field_number = copy.deepcopy(post_processesUnifield[i][
                0])
        
            else:

                valid_names = ""

                valid_set = list(fields_names.keys())

                for valid_name in valid_set:

                    valid_names += "'"+valid_name+"', "

                raise KeyError("'"+str(copy.deepcopy(
                post_processesUnifield[i][0]))+"' is not a valid name "+
                "for a field in the setting of post-processes in this "+
                "problem. Check out the valid names: "+valid_names)

        # Transforms the component dictionary into a live dictionary wi-
        # th proper methods and stuff

        post_processesUnifield[i] = post_processingSelectionSingleField(
        post_processesUnifield[i][1], context_class)

        # Gets the names of the processes for this field

        processes_names = list(post_processesUnifield[i].keys())

        # Adds the number of the field as key

        post_processesUnifield[i]["field number"] = field_number

        # Adds the pairs of field number and process names

        post_processesNamesList.append([i])

        for name in processes_names:

            post_processesNamesList[-1].append([field_number, name])

    # Returns the new live list of dictionaries

    return post_processesUnifield, post_processesNamesList