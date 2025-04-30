# Routine to store some methods to post-process solution inside the 
# pseudotime stepping methods

from dolfin import *

import source.post_processes.post_processes_functions as post_functions

########################################################################
#                   Post-process's classes templates                   #
########################################################################

# Defines a class to hold the data provided by the code as a context

class PostProcessContext:

    # Sets the common information provided to the class by the system

    def __init__(self, mesh, constitutive_model, dx, 
    domain_physGroupsNamesToTags):

        self.mesh = mesh

        self.constitutive_model = constitutive_model
        
        self.dx = dx

        self.physical_groupsList = set(dx.subdomain_data().array())
        
        self.domain_physGroupsNamesToTags = domain_physGroupsNamesToTags

# Defines a template for the post-processes' classes

class PostProcessMethod:

    def __init__(self, initialization_function, update_function,
    additional_information, code_providedInfo):

        # Each class must have four variables: the initialization func-
        # tion, i.e. the function that initializes all concerning varia-
        # bles to the post-process method; the update function, which 
        # runs the process proper; a list of required additional infor-
        # mation names, which will be used to retrieve information given
        # in dictionaries by the user, thus, the names are their keys; 
        # a list of variables given by the code, that are retrieven from
        # the context class

        self.initialization_function = initialization_function

        self.update_function = update_function

        self.additional_information = additional_information

        self.code_providedInfo = code_providedInfo

########################################################################
#                      Post-processing tools list                      #
########################################################################

# Sets a class for the method to save the field

class SaveField(PostProcessMethod):

    def __init__(self, context: PostProcessContext):

        # Initializes the parent template class and already passes to it
        # the initialization function, the update function, the additio-
        # nal data names, and the code-provided information

        super().__init__(post_functions.initialize_fieldSaving, 
        post_functions.update_fieldSaving, ["directory path", 
        "file name"], [])

# Sets a class for the method to save the stress field

class SaveStressField(PostProcessMethod):

    def __init__(self, context: PostProcessContext):

        # Initializes the parent template class and already passes to it
        # the initialization function, the update function, the additio-
        # nal data names, and the code-provided information

        super().__init__(post_functions.initialize_cauchyStressSaving, 
        post_functions.update_cauchyStressSaving, ["directory path", 
        "file name", "polynomial degree"], [context.mesh, 
        context.constitutive_model, context.dx, 
        context.physical_groupsList, context.domain_physGroupsNamesToTags])

# Sets a class for the method to homogenize a field

class HomogenizeField(PostProcessMethod):

    def __init__(self, context: PostProcessContext):

        # Initializes the parent template class and already passes to it
        # the initialization function, the update function, the additio-
        # nal data names, and the code-provided information

        super().__init__(post_functions.initialize_fieldHomogenization, 
        post_functions.update_fieldHomogenization, ["directory path", 
        "file name", "subdomain"], [context.dx])

# Sets a class for the method to homogenize the gradient of a field

class HomogenizeFieldsGradient(PostProcessMethod):

    def __init__(self, context: PostProcessContext):

        # Initializes the parent template class and already passes to it
        # the initialization function, the update function, the additio-
        # nal data names, and the code-provided information

        super().__init__(post_functions.initialize_gradientFieldHomogenization, 
        post_functions.update_gradientFieldHomogenization, ["directory"+
        " path", "file name", "subdomain"], [context.dx])