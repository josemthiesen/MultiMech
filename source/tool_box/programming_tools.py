# Routine to store decorators and programming utilities

import functools

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