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
elements_dictionary, mesh_dataClass, bilinear_form=0.0, linear_form=0.0, 
boundary_conditions=None, fluctuation_field=False):
    
    # Verifies if the multiscale_BCsDict is a dictionary

    if not isinstance(multiscale_BCsDict, dict):

        raise TypeError("The dictionary of multiscale boundary conditi"+
        "ons must be a dictionary, whose keys are the constrained fiel"+
        "ds' names\nand the values are dictionaries of boundary condti"+
        "ons' names and macro information")
    
    # Tests whether it is a dictionary of dictionaries

    if not isinstance(multiscale_BCsDict[list(multiscale_BCsDict.keys())
    [0]], dict):
        
        # Gives a generic field name and encapsulate it as a field dic-
        # tionary

        multiscale_BCsDict = {"generic_field": multiscale_BCsDict}

    # Tests whether any of the boundary condition is in factr a periodic
    # boundary condition

    for field_name, BC_type in multiscale_BCsDict.items():

        try: 

            if BC_type["boundary condition"][0:8]=="Periodic":

                # Updates the fluctuation field flag to true, because 
                # the periodic boundary condition must be coded with tje
                # fluctuation field
                
                fluctuation_field = True   

                # Checks the boundary condition's order

                BC_order = ""

                if BC_type["boundary condition"][-9:-1]=="stOrderB":

                    BC_order = "FirstOrderBC"

                elif BC_type["boundary condition"][-9:-1]=="ndOrderB":

                    BC_order = "SecondOrderBC"

                else:

                    raise NameError("The homogenization order of the '"+
                    BC_type["boundary condition"]+"' boundary conditio"+
                    "n could not be resolved. It does not end with 'Fi"+
                    "rstOrderBC' nor 'SecondOrderBC'")

                # Converts all boundary conditions to periodic, because
                # one field to be periodic requires all of them to be 
                # too

                for field_name in multiscale_BCsDict.keys():

                    multiscale_BCsDict[field_name]["boundary condition"
                    ] = "Periodic"+BC_order

        except:

            raise NameError("The dictionary of multiscale boundary con"+
            "ditions has a boundary condition dictionary attached to t"+
            "he field '"+str(field_name)+"' that does not have the key"+
            " 'boundary condition' or the corresponding value does not"+
            " have at least 8 characters.")

    # Gets the names of the fields from the keys of the dictionary of e-
    # lements

    fields_names = [field_name for field_name in elements_dictionary.keys(
    )]

    # Initializes the list of macro quantities classes as an empty list

    macro_quantitiesClasses = []

    # Initializes the inverse of the volume

    volume_inverse = None

    # Initializes a dictionary of field corrections, due to having or 
    # not the BVP been defined using a fluctuation field instead of the 
    # complete field

    field_corrections = dict()

    # Iterates through the fields

    for field_name, field_BC in multiscale_BCsDict.items():

        # Checks if this field has as value a dictionary with a boundary
        # condition key and a macro information key

        if not ("boundary condition" in field_BC):

            raise KeyError("The dictionary of multiscale boundary cond"+
            "itions has the fields' names as keys and another dictiona"+
            "ry as value.\nThis inner dictionary must have the key 'bo"+
            "undary condition'")

        if not ("macro information" in field_BC):

            raise KeyError("The dictionary of multiscale boundary cond"+
            "itions has the fields' names as keys and another dictiona"+
            "ry as value.\nThis inner dictionary must have the key 'ma"+
            "cro information'")

        # Constructs a dictionary of classes of multiscale boundary con-
        # ditions that are implemented. If you implement any new bounda-
        # ry condition, the code will automatically retrieve it using 
        # the inspect functionality if the class you created is a child 
        # of the class BCsClassTemplate.
        # Retrieves at last just the boundary condition because the dis-
        # patch functions returns a dictionary of methods

        multiscale_BCsDict[field_name] = programming_tools.dispatch_processes(
        field_BC["boundary condition"], multiscale_classes, 
        multiscale_classes.BCsClassTemplate, reserved_classes=[
        multiscale_classes.BCsClassTemplate], class_input=(field_name,
        fields_names, elements_dictionary, mesh_dataClass, field_BC["m"+
        "acro information"], macro_quantitiesClasses, volume_inverse,
        fluctuation_field))[field_BC["boundary condition"]]

        # Recovers the elements dictionary, the fields names, the macro
        # quantities list of classes, and the inverse of the volume

        elements_dictionary = multiscale_BCsDict[field_name
        ].elements_dictionary

        fields_names = multiscale_BCsDict[field_name].fields_names

        macro_quantitiesClasses = multiscale_BCsDict[field_name
        ].macro_quantitiesClasses

        volume_inverse = multiscale_BCsDict[field_name].volume_inverse

        # Recovers the correction of the primal field. If the BVP is 
        # constructed using the fluctuation field, a correction (linear
        # or quadratic) must be added. Constructs the function space to
        # project the solution with the correction later for visualiza-
        # tion and post-processing

        if fluctuation_field:

            field_corrections[field_name] = multiscale_BCsDict[
            field_name].field_correction

    # Constructs the mixed element using the order of fields from the 
    # list of field names, then, creates the monolithic function space

    mixed_element = MixedElement([elements_dictionary[field_name] for (
    field_name) in (fields_names)])

    monolithic_functionSpace = FunctionSpace(mesh_dataClass.mesh, 
    mixed_element)

    # Creates the function for the updated solution, i.e. the vector of 
    # parameters. Then, splits into the individual fields. As the mono-
    # lithic function space has an specific order of fields, retrieves 
    # the individual fields following this order. Does the same for the
    # variations

    trial_functions = TrialFunction(monolithic_functionSpace)

    monolithic_solution = Function(monolithic_functionSpace)

    solution_functions = split(monolithic_solution)

    variation_functions = split(TestFunction(monolithic_functionSpace))

    # Organize them into dictionaries

    solution_fields = dict()

    variation_fields = dict()

    fields_namesDict = dict()

    for i in range(len(fields_names)):

        field_name = fields_names[i]

        # Retrives the fields and adds the corrections due to using the
        # fluctuation field or not

        if fluctuation_field and (field_name in field_corrections):

            solution_fields[field_name] = (solution_functions[i]
            +field_corrections[field_name][0])

        else:

            solution_fields[field_name] = solution_functions[i]

        variation_fields[field_name] = variation_functions[i]

        fields_namesDict[field_name] = i

    # With the constructed fields, both trial and test functions, it is
    # time to build the variational forms and the boundary conditions

    # Initialize a counter of BCs classes

    counter_BCs = 0

    for field_BC in multiscale_BCsDict.values():

        # Uses the standard update method that is obligatory for all BCs
        # classes

        bilinear_form, linear_form, boundary_conditions = field_BC.update(
        bilinear_form, linear_form, boundary_conditions, solution_fields,
        variation_fields, counter_BCs, mesh_dataClass,
        monolithic_functionSpace)

        # Updates the counter of BCs

        counter_BCs += 1

    # Returns the variational forms, the boundary conditions, and the 
    # list of macro quantities classes

    return (bilinear_form, linear_form, boundary_conditions, 
    macro_quantitiesClasses, fields_namesDict, solution_fields, 
    variation_fields, trial_functions, monolithic_solution, 
    mixed_element, volume_inverse, field_corrections)