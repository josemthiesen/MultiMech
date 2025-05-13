# Routine to store methods to select and manage multiscale boundary con-
# ditions

from dolfin import *

import source.multiscale.multiscale_classes as multiscale_classes

import source.tool_box.programming_tools as programming_tools

# Defines a function to select and apply the multiscale boundary condi-
# tions. The boundary conditions are given by dictionaries of dictiona-
# ries: the outer dictionary is of fields, i.e. there is a key for the 
# name of each field, whereas the inner dictionary is a dictionary of 
# multiscale boundary conditions classes, the keys must match the names
# of the available classes of boundary conditions in multiscale_classes

@programming_tools.optional_argumentsInitializer({'boundary_conditions': 
lambda: []})

def select_multiscaleBoundaryConditions(multiscale_BCsDict,
elements_dictionary, mesh_dataClass, bilinear_form, linear_form, 
boundary_conditions=None):

    # Gets the names of the fields from the keys of the dictionary of e-
    # lements

    fields_names = [field_name for field_name in elements_dicionary.keys(
    )]

    # Initializes the list of macro quantities classes as an empty list

    macro_quantitiesClasses = []

    # Iterates through the fields

    for field_name, field_BCs in multiscale_BCsDict.items():

        # Constructs a dictionary of classes of multiscale boundary con-
        # ditions that are implemented. If you implement any new bounda-
        # ry condition, the code will automatically retrieve it using 
        # the inspect functionality if the class you created is a child 
        # of the class BCsClassTemplate

        multiscale_BCsDict[field_name] = programming_tools.dispatch_processes(
        field_BCs, multiscale_classes, 
        multiscale_classes.BCsClassTemplate, reserved_classes=[
        multiscale_classes.BCsClassTemplate], class_input=(field_name,
        field_name+" gradient", fields_names, elements_dictionary, 
        mesh_dataClass, macro_quantitiesClasses))

        # Recovers the elements dictionary, the fields names, and the 
        # macro quantities list of classes

        elements_dictionary = multiscale_BCsDict[field_name
        ].elements_dictionary

        fields_names = multiscale_BCsDict[field_name].fields_names

        macro_quantitiesClasses = multiscale_BCsDict[field_name
        ].macro_quantitiesClasses

    # Constructs the mixed element using the order of fields from the 
    # list of field names, then, creates the monolithic function space

    mixed_element = [elements_dictionary[field_name] for field_name in (
    fields_names)]

    monolithic_functionSpace = FunctionSpace(mesh_dataClass.mesh, 
    mixed_element)

    # Defines the trial functions

    delta_solution = TrialFunction(monolithic_functionSpace) 

    # Creates the function for the updated solution, i.e. the vector of 
    # parameters. Then, splits into the individual fields. As the mono-
    # lithic function space has an specific order of fields, retrieves 
    # the individual fields following this order. Does the same for the
    # variations

    solution_functions = split(Function(monolithic_functionSpace))

    variation_functions = split(TestFunction(monolithic_functionSpace))

    # Organize them into dictionaries

    solution_fields = dict()

    variation_fields = dict()

    fields_namesDict = dict()

    for i in range(len(fields_names)):

        solution_fields[fields_names[i]] = solution_functions[i]

        variation_fields[fields_names[i]] = variation_functions[i]

        fields_namesDict[fields_names[i]] = i

    # With the constructed fields, both trial and test functions, it is
    # time to build the variational forms and the boundary conditions

    # Initialize a counter of BCs classes

    counter_BCs = 0

    for field_BCs in multiscale_BCsDict.values():

        # Uses the standard update method that is obligatory for all BCs
        # classes

        bilinear_form, linear_form, boundary_conditions = fields_BCs.update(
        bilinear_form, linear_form, boundary_conditions, solution_fields,
        variation_fields, counter_BCs, mesh_dataClass)

        # Updates the counter of BCs

        counter_BCs += 1

    # Returns the variational forms, the boundary conditions, and the 
    # list of macro quantities classes

    return (bilinear_form, linear_form, boundary_conditions, 
    macro_quantitiesClasses, fields_namesDict, solution_fields, 
    variation_fields)