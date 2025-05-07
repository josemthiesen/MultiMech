# Routine to store some methods to post-process solution inside the 
# pseudotime stepping methods

import copy

import inspect

import source.post_processes.post_processes_classes as post_classes

########################################################################
#             Post-processes' classes automatic retriever              #
########################################################################

# Defines a function to get the classes of the post_classes module

def get_allProcesses():

    # Initializes a dictionary of classes for the post-process methods

    methods_classesDict = dict()

    # Checks for the classes in the post_classes module using inspect 
    # library

    for class_name, class_object in inspect.getmembers(post_classes,
    inspect.isclass):
        
        # Checks if the class is a subclass of the PostProcessMethod pa-
        # rent class

        if (issubclass(class_object, post_classes.PostProcessMethod) and (
        class_object is not post_classes.PostProcessMethod) and (
        class_object is not post_classes.PostProcessContext)):
            
            # Appends this subclass

            methods_classesDict[class_name] = class_object

    # Returns the lis of classes

    return methods_classesDict

# Defines a function to initializes the needed methods' classes

def dispatch_processes(methods_names, context_class):

    # Initializes a dictionary of classes that will be used in the post-
    # processes

    post_processClasses = dict()

    # Gets the methods classes of this file as a dictionary, the keys a-
    # re the names of the classes whereas the values are the classes ob-
    # jects

    methods_classesDict = get_allProcesses()

    available_methodsNames = list(methods_classesDict.keys())

    # Iterates through the methods names

    for method_name in methods_names:

        # If the name is in the list of available methods, updates the
        # dictionary of classes used in the post-processing step

        if method_name in available_methodsNames:

            post_processClasses[method_name] = methods_classesDict[
            method_name](context_class)

        else:

            available_list = "\n"

            for name in available_methodsNames:

                available_list += "'"+name+"'\n"

            raise NameError("'"+str(method_name)+"' is not an availabl"+
            "e post-processing method. Find one in the list:"+
            available_list)
        
    # Returns the dictionary of classes that will be used in the post-
    # processing step

    return post_processClasses

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

    available_processes = dispatch_processes(post_processes, 
    context_class)

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
# tion. The post_processesMultifield is a dictionary of processes that
# must receive all or a set of fields as input

def post_processingSelectionMultipleFields(post_processesUnifield,
n_fields, context_class):
    
    # Iterates through the number of fields

    for i in range(len(post_processesUnifield)):

        # Transforms the component dictionary into a live dictionary wi-
        # th proper methods and stuff

        post_processesUnifield[i] = post_processingSelectionSingleField(
        post_processesUnifield[i], context_class)

    # Returns the new live list of dictionaries

    return post_processesUnifield