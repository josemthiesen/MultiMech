# Routine to store decorators and programming utilities

import functools

import inspect

########################################################################
#                         Return-handling tools                        #
########################################################################

# Defines a function to get a particular variable from the result of a 
# function

def get_result(whole_result, variable_name):

    # Verifies whether the whole result is a dictionary

    if isinstance(whole_result, dict):

        # Tries to read the variable name as a key

        try:

            return whole_result[variable_name]
        
        except:

            raise KeyError("The name '"+str(variable_name)+"' is not a"+
            " key of the result dictionary of this function. Do you me"+
            "an one of the following possible keys? "+str(set(
            whole_result.keys()))[1:-1])
        
    # Verifies if it is a tuple

    if isinstance(whole_result, tuple):

        # Verifies if it is a named tuple
        
        if hasattr(whole_result, '_fields'):
        
            try:

                return getattr(whole_result, variable_name)
            
            except:

                raise AttributeError("The named tuple returned by the "+
                "function does not have a variable named "+str(
                variable_name)+". This named tuple has the following a"+
                "ttributes:\n"+str(whole_result._fields)[1:-1])
            
        # If it is not a named tuple, raise an exception

        else:

            raise TypeError("The result of the function is a tuple but"+
            " it is not a named tuple. Thus, the variable name '"+str(
            variable_name)+"' cannot be used")
        
    # Simply gives the result
        
    else:

        return whole_result

########################################################################
#                         Memory-handling tools                        #
########################################################################

# Defines a function to handle optional mutable arguments, to avoid re-
# writing in memory during multiple method calls. Use it as a decorator:
#
# @optional_argumentsInitializer({argument1: value1, argument2: value2})
# def function_example(arguments, argument1=None, argument2=None):
#     ...
# The value must be function so that it gives fresh values at every time

def optional_argumentsInitializer(default_argumentsDictionary):

    # Checks if all the values of the dictionary are callable, i.e. if
    # they are all functions

    for key, value in default_argumentsDictionary.items():

        if not callable(value):

            raise TypeError("All values in the dictionary of the decor"+
            "ator to get the default values for optional arguments mus"+
            "t be callable. Use lambda functions or functions.")

    # Defines the decorator function. It receives the decorated function
    # as argument

    def decorator(decorated_function):

        @functools.wraps(decorated_function)

        # Defines the wrapper to get the arguments and check if any of 
        # them is in the list of arguments to be initialized as prescri-
        # bed

        def wrapper(*positional_arguments, **keyword_argument):

            # Gets the names of the arguments using internal Python syn-
            # tax. Then, turns the arguments as a list

            arguments_names = decorated_function.__code__.co_varnames

            positional_arguments = list(positional_arguments)

            # Iterates through the arguments names and gets the index of
            # each name

            for index, name in enumerate(arguments_names):

                # Checks if this argument is in the names of arguments
                # that are prescribed by default

                if name in default_argumentsDictionary:

                    # If the index is still less than the number of po-
                    # sitionally given arguments

                    if index<len(positional_arguments):

                        # Checks if the given value is None; if so, 
                        # substitutes it by the default value

                        if positional_arguments[index] is None:

                            positional_arguments[index] = (
                            default_argumentsDictionary[name]())

                    # If it is a keyword argument

                    else:

                        # Checks if the given value is None; if so, 
                        # substitutes it by the default value

                        if keyword_argument.get(name, None) is None:

                            keyword_argument[name] = (
                            default_argumentsDictionary[name]())

            return decorated_function(*positional_arguments, 
            **keyword_argument)
        
        return wrapper
    
    return decorator

########################################################################
#                Classes from file automatic retriever                 #
########################################################################

# Defines a function to get the classes of a python module (
# searched_file) which are subclasses of a parent class. The 
# reserved_classes is a list of classes that shall not be given as a 
# process for actual evaluation, they can be parent classes for example

@optional_argumentsInitializer({'reserved_classes': lambda: []})

def get_allProcesses(searched_file, parent_class, reserved_classes=None):

    # Initializes a dictionary of classes for the classes found in a fi-
    # le

    methods_classesDict = dict()

    # Checks for the classes in the module using inspect library

    for class_name, class_object in inspect.getmembers(searched_file,
    inspect.isclass):
        
        # Checks if the class is a subclass of the PostProcessMethod pa-
        # rent class

        if (issubclass(class_object, parent_class) and (not (class_object
        ) in reserved_classes)):
            
            # Appends this subclass

            methods_classesDict[class_name] = class_object

    # Returns the lis of classes

    return methods_classesDict

# Defines a function to initializes the needed methods' classes

@optional_argumentsInitializer({'reserved_classes': lambda: []})

def dispatch_processes(methods_names, searched_file, parent_class, 
class_input=None, reserved_classes=None):

    # Assures methods_names is a list

    if not (isinstance(methods_names, list) or isinstance(methods_names, 
    tuple)):

        methods_names = [methods_names]

    # Initializes a dictionary of classes that will be used in the post-
    # processes

    methods_classes = dict()

    # Gets the methods classes of this file as a dictionary, the keys a-
    # re the names of the classes whereas the values are the classes ob-
    # jects

    methods_classesDict = get_allProcesses(searched_file, parent_class,
    reserved_classes=reserved_classes)

    available_methodsNames = list(methods_classesDict.keys())

    # Iterates through the methods names

    if class_input is None:

        for method_name in methods_names:

            # If the name is in the list of available methods, updates 
            # the dictionary of classes

            if method_name in available_methodsNames:

                methods_classes[method_name] = methods_classesDict[
                method_name]()

            else:

                available_list = "\n"

                for name in available_methodsNames:

                    available_list += "'"+name+"'\n"

                raise NameError("'"+str(method_name)+"' is not an avai"+
                "lable class method. Find one in the list:"+
                available_list)

    elif isinstance(class_input, tuple):

        for method_name in methods_names:

            # If the name is in the list of available methods, updates 
            # the dictionary of classes

            if method_name in available_methodsNames:

                methods_classes[method_name] = methods_classesDict[
                method_name](*class_input)

            else:

                available_list = "\n"

                for name in available_methodsNames:

                    available_list += "'"+name+"'\n"

                raise NameError("'"+str(method_name)+"' is not an avai"+
                "lable class method. Find one in the list:"+
                available_list)

    else:

        for method_name in methods_names:

            # If the name is in the list of available methods, updates 
            # the dictionary of classes

            if method_name in available_methodsNames:

                methods_classes[method_name] = methods_classesDict[
                method_name](class_input)

            else:

                available_list = "\n"

                for name in available_methodsNames:

                    available_list += "'"+name+"'\n"

                raise NameError("'"+str(method_name)+"' is not an avai"+
                "lable class method. Find one in the list:"+
                available_list)
        
    # Returns the dictionary of classes that will be used given the na-
    # mes and the available classes

    return methods_classes