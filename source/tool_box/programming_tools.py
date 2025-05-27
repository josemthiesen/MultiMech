# Routine to store decorators and programming utilities

import functools

import inspect

from collections import OrderedDict

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
#                           Dictionary tools                           #
########################################################################

# Defines a function to swap the keys of a dictionary. Uses a list of 
# lists to change some or all keys; each sublist is a pair of old and 
# new keys

def change_dictionaryKeys(original_dictionary, old_keysList):

    for old_key, new_key in old_keysList:

        # Tests if the key is in the dictionary

        if old_key in original_dictionary:

            # Adds the new key

            original_dictionary[new_key] = original_dictionary[old_key]

            # Deletes the old key

            original_dictionary.pop(old_key)

        # If this key is not found, raises an error message

        else:

            raise KeyError("The key '"+str(old_key)+"' cannot be delet"+
            "ed and swapped by the new key '"+str(new_key)+"' because "+
            "the old keys does not exist in the dictionary.")

    return original_dictionary

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
#                Classes from file automatic retriever                 #
########################################################################

# Defines a function to get the entities of a python module (
# searched_file) which are subclasses of a parent class. The 
# reserved_entities is a list of entities that shall not be given as a 
# process for actual evaluation, they can be parent classes for example

@optional_argumentsInitializer({'reserved_entities': lambda: []})

def get_allProcesses(searched_file, parent_class=None, reserved_entities=
None, searched_entity="classes"):
    
    # Sets the stuff you want to search for

    if searched_entity=="classes":

        searched_entity = inspect.isclass

    elif searched_entity=="functions":

        searched_entity = inspect.isfunction

    else:

        raise NameError("The entitities '"+str(searched_entity)+"' can"+
        "not be searched. You can select from either 'classes' or 'fun"+
        "ctions'")

    # Initializes a dictionary of entities for the entities found in a 
    # file

    entities_dict = dict()

    # If the searched entity is a function or the searched classes are 
    # not children of some parent class

    if searched_entity==inspect.isfunction or (parent_class is None):

        # Checks for the entities in the module using inspect library

        for entity_name, entity_object in inspect.getmembers(
        searched_file, searched_entity):
            
            # Checks if the class is a subclass of the parent class

            if not (entity_object in reserved_entities):
                
                # Appends this subclass

                entities_dict[entity_name] = entity_object

    # If the searched entity is a class and the parent class is not None

    elif searched_entity==inspect.isclass:

        # Checks for the entities in the module using inspect library

        for entity_name, entity_object in inspect.getmembers(
        searched_file, searched_entity):
            
            # Checks if the class is a subclass of the parent class

            if (issubclass(entity_object, parent_class) and (not (
            entity_object in reserved_entities))):
                
                # Appends this subclass

                entities_dict[entity_name] = entity_object

    # Returns the lis of classes

    return entities_dict

# Defines a function to initializes the needed methods' classes

@optional_argumentsInitializer({'reserved_entities': lambda: []})

def dispatch_classes(methods_names, searched_file, parent_class=None, 
class_input=None, reserved_entities=None):

    # Assures methods_names is a list, or a tuple or a dictionary

    if not (isinstance(methods_names, list) or isinstance(methods_names, 
    tuple) or isinstance(methods_names, dict)):

        methods_names = [methods_names]

    # Initializes a dictionary of classes that will be used later

    methods_classes = dict()

    # Gets the methods classes of this file as a dictionary, the keys a-
    # re the names of the classes whereas the values are the classes ob-
    # jects

    methods_classesDict = get_allProcesses(searched_file, parent_class=
    parent_class, reserved_entities=reserved_entities)

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

# Defines a function to initializes the needed methods' functions

@optional_argumentsInitializer({'reserved_entities': lambda: []})

def dispatch_functions(methods_names, searched_file, reserved_entities=
None, fixed_inputVariablesDict=None, methods_functionsDict=None, 
return_list=False, return_singleFunction=False, all_argumentsFixed=False):

    # Assures methods_names is a list, or a tuple or a dictionary

    if not (isinstance(methods_names, list) or isinstance(methods_names, 
    tuple) or isinstance(methods_names, dict)):

        methods_names = [methods_names]

    # Initializes a dictionary of functions that will be used later

    methods_functions = OrderedDict()

    # Gets the functions of this file as a dictionary, the keys are the
    # names of the functions whereas the values are the functions objects

    if methods_functionsDict is None:

        methods_functionsDict = get_allProcesses(searched_file, 
        reserved_entities=reserved_entities, searched_entity="functions")

    available_methodsNames = list(methods_functionsDict.keys())

    # Iterates through the methods names

    if fixed_inputVariablesDict is None:

        # Iterates through the methods

        for method_name in methods_names:

            # If the name is in the list of available methods, updates 
            # the dictionary of classes

            if method_name in available_methodsNames:

                methods_functions[method_name] = methods_functionsDict[
                method_name]

            else:

                available_list = "\n"

                for name in available_methodsNames:

                    available_list += "'"+name+"'\n"

                raise NameError("'"+str(method_name)+"' is not an avai"+
                "lable function method. Find one in the list:"+
                available_list)

    else:

        # Verifies if the dictionary of input variables to be fixed is 
        # indeed a dictionary

        if not isinstance(fixed_inputVariablesDict, dict):

            raise TypeError("The dictionary of fixed input variables m"+
            "ust be a dictionary")

        # Iterates through the methods

        for method_name in methods_names:

            # Verifies if the dictionary of input variables to be fixed
            # has this method as a key

            fixed_arguments = []

            if not (method_name in available_methodsNames):

                raise KeyError("The '"+str(method_name)+"' is not an a"+
                "vailable method. Check out the available methods: "+str(
                available_methodsNames))

            if isinstance(fixed_inputVariablesDict, dict):

                # Checks if the dictionary of fixed input variables is a
                # dictionary of dictionaries

                if isinstance(fixed_inputVariablesDict[list(
                fixed_inputVariablesDict.keys())[0]], dict):

                    if not (method_name in fixed_inputVariablesDict):

                        raise KeyError("The dictionary of fixed input vari"+
                        "ables does not have the method '"+str(method_name)+
                        "' as a key. The available methods are: "+str(
                        fixed_inputVariablesDict.keys()))

                    fixed_arguments = fixed_inputVariablesDict[method_name]
                
                else:

                    fixed_arguments = fixed_inputVariablesDict
            
            else:

                raise TypeError("The corresponding value for the metho"+
                "d '"+str(method_name)+"' in the dictionary of fixed i"+
                "nput variables is not a dictionary")

            # If the name is in the list of available methods, updates 
            # the dictionary of classes

            if method_name in available_methodsNames:

                # If there is no fixed arguments, just gets the function

                if len(fixed_arguments.keys())==0:

                    methods_functions[method_name] = methods_functionsDict[
                    method_name]

                # Otherwise, separates the fixed arguments from the free
                # ones

                else:

                    # Sets the new function as a driver, i.e. fixing so-
                    # me arguments given by the fixed arguments dictio-
                    # nary
                    
                    methods_functions[method_name] = driver_maker(
                    methods_functionsDict[method_name], method_name, 
                    fixed_arguments, all_argumentsFixed=
                    all_argumentsFixed)

            else:

                available_list = "\n"

                for name in available_methodsNames:

                    available_list += "'"+name+"'\n"

                raise NameError("'"+str(method_name)+"' is not an avai"+
                "lable function method. Find one in the list:"+
                available_list)
            
    # Verifies if the output must be a list instead of a dictionary. When
    # usage is direct

    if return_list:

        methods_functions = list(methods_functions.values())

        if return_singleFunction and len(methods_functions)==1:

            methods_functions = methods_functions[0]

    return methods_functions, methods_functionsDict

# Defines a maker of driver functions by fixing the arguments of some 
# other function

def driver_maker(method_function, method_name, fixed_arguments, 
all_argumentsFixed=False):

    # Gets the names of the arguments of the function that is being as-
    # signed

    arguments_names = inspect.signature(method_function).parameters

    # Separates the positional arguments from the keyword ones

    positional_arguments = [name for name, param in (
    arguments_names.items()) if param.default is inspect.Parameter.empty]

    keyword_arguments = [name for name, param in arguments_names.items(
    ) if param.default is not inspect.Parameter.empty]

    # Initializes a dictionary of fixed arguments that are really argu-
    # ments of the function

    functions_fixedArguments = dict()

    # If more fixed arguments were provided, thus, enabling all the ar-
    # guments to be picked from a larger pool of arguments

    if all_argumentsFixed:

        # Verifies if any of the positional arguments was not provided 
        # in the fixed arguments

        for fixed_argument in positional_arguments:

            if not (fixed_argument in fixed_arguments):

                raise NameError("The argument '"+str(fixed_argument)+
                "' was not provided in the list of fixed arguments. Ch"+
                "eck the list of fixed arguments: "+str(
                fixed_arguments.keys()))
    
            # If this argument's name is in the general dictionary of 
            # fixed arguments, updates the functions own dictionary of
            # this kind
            
            else:
                
                functions_fixedArguments[fixed_argument] = (
                fixed_arguments[fixed_argument])

    else:

        # Verifies if any of the fixed arguments is not really an argu-
        # ment of the required function

        for fixed_argument in fixed_arguments:

            if not (fixed_argument in arguments_names):

                raise NameError("The fixed argument '"+str(
                fixed_argument)+"' is not an argument of the function "+
                "'"+str(method_name)+"'. Did you mean any of the follo"+
                "wing arguments? "+str(arguments_names))
    
            # If this argument's name is in the general dictionary of 
            # fixed arguments, updates the functions own dictionary of
            # this kind
            
            else:
                
                functions_fixedArguments[fixed_argument] = (
                fixed_arguments[fixed_argument])

    # Gets the free arguments from the difference with the list of fixed 
    # arguments. Again differentiates the positional from the keyword 
    # arguments

    free_positionalNames = [free_name for free_name in (
    positional_arguments) if free_name not in functions_fixedArguments]

    free_keywordNames = [free_name for free_name in (keyword_arguments
    ) if free_name not in functions_fixedArguments]

    # Creates the driver of the function using a wrapper to keep the o-
    # riginal function's information

    @functools.wraps(method_function)

    def driver_function(*arguments, **key_arguments):

        # Checks if the arguments given have the same length of the re-
        # quired arguments

        if len(arguments)!=len(free_positionalNames):

            raise ValueError("You've provided "+str(len(arguments))+" "+
            "positional arguments to the function '"+str(method_name)+
            "' whereas only "+str(len(free_positionalNames))+" are req"+
            "uired. To check, read the list of positional arguments ah"+
            "ead: "+str(free_positionalNames))

        # Checks if the key arguments are really keys

        for key in key_arguments.keys():

            if not (key in free_keywordNames):

                raise KeyError("The key '"+str(key)+"' is not keyword "+
                "argument of the function '"+str(method_name)+"'. The "+
                "following is a list of the keyword arguments of this "+
                "method: "+str(free_keywordNames))

        return method_function(**{free_positionalNames[i]: arguments[i
        ] for i in range(len(free_positionalNames))}, **key_arguments, 
        **functions_fixedArguments)
    
    return driver_function