# Routine to store variational formulation methods

from dolfin import *

import ufl_legacy

import copy

import source.tool_box.tensor_tools as tensor_tools

import source.tool_box.programming_tools as programming_tools

import source.tool_box.surface_loading_tools as surface_loading_tools

########################################################################
#                            Internal work                             #
########################################################################

# Defines a function to construct the variational form of a non-dissipa-
# tive solid mechanics problem. The constitutive model can be either a
# class (when the whole domain has only one constitutive model); or it 
# can be a dictionary, when the domain is heterogeneous. The keys of the
# dictionaries are the volumetric physical groups, whereas the values 
# are the constitutive model classes. This internal work is calculated 
# using the first Piola-Kirchhoff stress tensor

def hyperelastic_internalWorkFirstPiola(trial_function, test_function, 
constitutive_modelDictionary, mesh_dataClass):
    
    # Gets the physical groups from the domain mesh function

    physical_groupsList = set(mesh_dataClass.dx.subdomain_data().array())

    # Initializes the variational form of the inner work

    inner_work = 0.0

    # If the constitutive model is a dictionary, the domain is heteroge-
    # neous

    if isinstance(constitutive_modelDictionary, dict):

        # Iterates through the dictionary

        for physical_group, constitutive_model in (
        constitutive_modelDictionary.items()):
            
            # Verifies physical group consistency, i.e. if it exists and
            # if it is an integer or a tuple of integers

            physical_group = verify_physicalGroups(physical_group, 
            physical_groupsList, physical_groupsNamesToTags=
            mesh_dataClass.domain_physicalGroupsNameToTag)

            # Initializes objects for the stresses at the reference 
            # configuration

            first_piola = programming_tools.get_result(
            constitutive_model.first_piolaStress(trial_function), "fir"+
            "st_piola_kirchhoff")

            # If the physical group is indeed a list of tags, i.e. mul-
            # tiple regions are integrated using the same constitutive
            # model

            if isinstance(physical_group, list):

                # Iterates through the subdomains

                for sub_physicalGroup in physical_group:

                    print("The volume of the", sub_physicalGroup, "sub"+
                    "domain is:", assemble(1*mesh_dataClass.dx(
                    sub_physicalGroup)), "\n")

                    # Constructs the variational forms for the inner 
                    # work

                    inner_work += (inner(first_piola, grad(test_function
                    ))*mesh_dataClass.dx(sub_physicalGroup))

            else:

                print("The volume of the", physical_group, "subdomain "+
                "is:", assemble(1*mesh_dataClass.dx(physical_group)), 
                "\n")

                # Constructs the variational forms for the inner work

                inner_work += (inner(first_piola, grad(test_function))*
                mesh_dataClass.dx(physical_group))

    # If the constitutive model is not a dictionary, the domain is homo-
    # geneous

    else:

        print("The volume of the domain is:", assemble(1*
        mesh_dataClass.dx), "\n")

        # Initializes objects for the stresses at the reference configu-
        # ration

        first_piola = programming_tools.get_result(
        constitutive_modelDictionary.first_piolaStress(trial_function), 
        "first_piola_kirchhoff")

        # Constructs the variational forms for the inner work

        inner_work = (inner(first_piola, grad(test_function))*
        mesh_dataClass.dx)

    # Returns the inner work variational form

    if mesh_dataClass.verbose:

        print("Finishes creating the variational form of the inner wor"+
        "k done by the\nfirst Piola stress tensor in a Cauchy continuu"+
        "m medium\n")

    return inner_work

# Defines a function to construct the variational form of a non-dissipa-
# tive micropolar continuum problem. The constitutive model can be ei-
# ther a class (when the whole domain has only one constitutive model); 
# or it can be a dictionary, when the domain is heterogeneous. The keys 
# of the dictionaries are the volumetric physical groups, whereas the 
# values are the constitutive model classes. This internal work is cal-
# culated using the first Piola-Kirchhoff stress tensor and its couple
# stress

def hyperelastic_micropolarInternalWorkFirstPiola(
displacement_trialFunction, microrotation_trialFunction, 
displacement_testFunction, microrotation_testFunction,
constitutive_modelDictionary, mesh_dataClass):
    
    # Gets the physical groups from the domain mesh function

    physical_groupsList = set(mesh_dataClass.dx.subdomain_data().array(
    ))

    # Initializes the variational form of the inner work

    inner_work = 0.0

    # If the constitutive model is a dictionary, the domain is heteroge-
    # neous

    if isinstance(constitutive_modelDictionary, dict):

        # Iterates through the dictionary

        for physical_group, constitutive_model in (
        constitutive_modelDictionary.items()):
            
            # Verifies physical group consistency, i.e. if it exists and
            # if it is an integer or a tuple of integers

            physical_group = verify_physicalGroups(physical_group, 
            physical_groupsList, physical_groupsNamesToTags=
            mesh_dataClass.domain_physicalGroupsNameToTag)

            # Initializes objects for the stresses at the reference 
            # configuration

            result = constitutive_model.first_piolaStress([
            displacement_trialFunction, microrotation_trialFunction])

            first_piola = programming_tools.get_result(result, "first_"+
            "piola_kirchhoff")

            couple_firstPiola = programming_tools.get_result(result, 
            "couple_first_piola_kirchhoff")

            kirchhoff = programming_tools.get_result(
            constitutive_model.kirchhoff_stress([
            displacement_trialFunction, microrotation_trialFunction]),
            "kirchhoff")

            # If the physical group is indeed a list of tags, i.e. mul-
            # tiple regions are integrated using the same constitutive
            # model

            if isinstance(physical_group, list):

                # Iterates through the individual physical groups

                for sub_physicalGroup in physical_group:

                    # Constructs the variational forms for the inner 
                    # work of the first Piola-Kirchhoff stress

                    inner_work += (inner(first_piola, grad(
                    displacement_testFunction))*mesh_dataClass.dx(
                    sub_physicalGroup))

                    # Adds the parcel of the couple stress

                    inner_work += ((inner(couple_firstPiola, grad(
                    microrotation_testFunction))*mesh_dataClass.dx(
                    sub_physicalGroup))-(inner(kirchhoff, 
                    tensor_tools.skew_2OrderTensor(
                    microrotation_testFunction))*mesh_dataClass.dx(
                    sub_physicalGroup)))

            else:

                # Constructs the variational forms for the inner work of
                # the first Piola-Kirchhoff stress

                inner_work += (inner(first_piola, grad(
                displacement_testFunction))*mesh_dataClass.dx(
                physical_group))

                # Adds the parcel of the couple stress

                inner_work += ((inner(couple_firstPiola, grad(
                microrotation_testFunction))*mesh_dataClass.dx(
                physical_group))-(inner(kirchhoff, 
                tensor_tools.skew_2OrderTensor(
                microrotation_testFunction))*mesh_dataClass.dx(
                physical_group)))

    # If the constitutive model is not a dictionary, the domain is homo-
    # geneous

    else:

        # Initializes objects for the stresses at the reference configu-
        # ration

        result = constitutive_model.first_piolaStress([
        displacement_trialFunction, microrotation_trialFunction])

        first_piola = programming_tools.get_result(result, "first_piol"+
        "a_kirchhoff")

        couple_firstPiola = programming_tools.get_result(result, "coup"+
        "le_first_piola_kirchhoff")

        kirchhoff = programming_tools.get_result(
        constitutive_model.kirchhoff_stress([displacement_trialFunction, 
        microrotation_trialFunction]), "kirchhoff")

        # Constructs the variational forms for the inner work of the
        # first Piola-Kirchhoff stress

        inner_work += (inner(first_piola, grad(displacement_testFunction
        ))*mesh_dataClass.dx)

        # Adds the parcel of the couple stress

        inner_work += ((inner(couple_firstPiola, grad(
        microrotation_testFunction))*mesh_dataClass.dx)-(inner(kirchhoff, 
        tensor_tools.skew_2OrderTensor(microrotation_testFunction))*
        mesh_dataClass.dx))

    # Returns the inner work variational form

    if mesh_dataClass.verbose:

        print("Finishes creating the variational form of the inner wor"+
        "k done by the\nfirst Piola stress tensor and its couple stres"+
        "s in a micropolar continu-\num medium\n")

    return inner_work

########################################################################
#               Work done by Neumann boundary conditions               #
########################################################################

# Defines a function to construct the variational form of the work done
# by the traction vector in the reference configuration given a dictio-
# nary of traction loads, where the keys are the corresponding boundary
# physical groups and the values are the traction loads

def traction_work(traction_dictionary, field, field_variation, 
mesh_dataClass, neumann_loads):

    # Gets the physical groups tags

    physical_groupsTags = set(mesh_dataClass.ds.subdomain_data().array(
    ))

    # Initializes the variational form

    traction_form = 0.0

    # Initializes a dictionary of load-generating functions from the 
    # surface_loading_tools file

    methods_functionsDict = None

    methods_functionsDict = programming_tools.dispatch_functions([], 
    surface_loading_tools, methods_functionsDict=methods_functionsDict)[1]

    # Initializes the dictionary of fixed arguments for the loading 
    # functions

    fixed_arguments = {"field": field, "mesh_dataClass": mesh_dataClass,
    "field_variation": field_variation}

    # Iterates through the dictionary

    for physical_group, traction in traction_dictionary.items():

        # Verifies if the traction is a list, to add multiple loads to a
        # single physical group

        if isinstance(traction, list):

            # Iterates through the loads

            for load in traction:

                # Updates the traction form

                traction_form, neumann_loads = set_tractionIntegration(
                traction_form, traction, physical_group, 
                physical_groupsTags, mesh_dataClass, fixed_arguments, 
                methods_functionsDict, field_variation, neumann_loads)

        # If the traction is not a list, updates the variational form 
        # directly

        else:
            
            traction_form, neumann_loads = set_tractionIntegration(
            traction_form, traction, physical_group, physical_groupsTags, 
            mesh_dataClass, fixed_arguments, methods_functionsDict, 
            field_variation, neumann_loads)

    # Returns the variational form

    if mesh_dataClass.verbose:

        print("Finishes creating the variational form of the work done"+
        " by the traction on the boundary\n")

    return traction_form, neumann_loads

# Defines a function to integrate the traction variational form

def set_tractionIntegration(traction_form, traction, physical_group,
physical_groupsTags, mesh_dataClass, fixed_arguments, 
methods_functionsDict, field_variation, neumann_loads):

    # Checks if this traction is a dictionary with instructions

    if isinstance(traction, dict):

        # Checks if there is a load case name

        if not ("load case" in traction):

            raise KeyError("There is no key 'load case' in the diction"+
            "ary of traction for the physical group '"+str(
            physical_group)+"'. This key must be in to signal which au"+
            "tomatically-generated load case must be used. The followi"+
            "ng keys have been found though: "+str(traction.keys()))

        # Assembles the input of arguments for the traction vector
        # building method. Adds first the values given in fixed_ar-
        # guments

        method_arguments = {key: value for key, value in (
        fixed_arguments.items())}

        for key, value in traction.items():

            if key!="load case":

                method_arguments[key] = value

        # Dispatches the function and calls it right away

        traction, neumann_load = programming_tools.dispatch_functions(
        traction["load case"], None, fixed_inputVariablesDict=
        method_arguments, methods_functionsDict=methods_functionsDict, 
        return_list=True, return_singleFunction=True, all_argumentsFixed=
        True)[0]()

        # Appends the neumann_load to the list of time controls

        neumann_loads.append(neumann_load)

    # Verifies if the traction is, then, a fenics format

    elif ((not isinstance(traction, Expression)) and (not isinstance(
    traction, Constant)) and (not isinstance(traction, 
    ufl_legacy.core.expr.Expr))):

        raise TypeError("The traction vector, if not defined as a dict"+
        "ionary of instructions to use built-in load cases, must be de"+
        "fined as a dolfin format, either a Constant, an Expression, o"+
        "r as_vector. In the physical group '"+str(physical_group)+"',"+
        " the traction provided was: "+str(traction))

    # Verifies if this physical group is indeed in ds
        
    physical_group = verify_physicalGroups(physical_group, 
    physical_groupsTags, physical_groupsNamesToTags=
    mesh_dataClass.boundary_physicalGroupsNameToTag)

    if isinstance(physical_group, list):

        for sub_physicalGroup in physical_group:

            print("The physical group "+str(sub_physicalGroup)+" h"+
            "as an area of "+str(assemble(1*mesh_dataClass.ds(
            sub_physicalGroup)))+"\n")

            traction_form += (dot(traction, field_variation)*
            mesh_dataClass.ds(sub_physicalGroup))

    else:

        print("The physical group "+str(physical_group)+" has an a"+
        "rea of "+str(assemble(1*mesh_dataClass.ds(physical_group)))
        +"\n")

        traction_form += (dot(traction, field_variation)*
        mesh_dataClass.ds(physical_group))

    return traction_form, neumann_loads

########################################################################
#                              Utilities                               #
########################################################################

# Defines a function to verify if a physical group is consistent and if
# it is the whole list of existing physical groups

@programming_tools.optional_argumentsInitializer({('physical_groupsNam'+
'esToTags'): lambda: dict()})

def verify_physicalGroups(physical_group, physical_groupsList, 
physical_groupsNamesToTags=None, throw_error=True):
            
    # If the key of the dictionary is a string, it is the physical
    # group's name. Hence, converts it to its corresponding number tag
    
    if isinstance(physical_group, str):

        try:

            print("Transforms the physical group name '"+physical_group+
            "' to the tag:")

            physical_group = physical_groupsNamesToTags[physical_group]

            print(physical_group, "\n")

        except:

            if throw_error:

                raise KeyError("The physical group name '"+physical_group
                +"' was used in a variational form, but it does not ex"+
                "ist in the dictionary of physical groups' names to ta"+
                "gs. This dictionary has the following keys and values"+
                ": "+str(physical_groupsNamesToTags))

            else:

                return physical_group

    # Tuples can be used as physical groups to integrate over multiple 
    # physical groups simultaneously

    elif (not isinstance(physical_group, int)) and (not isinstance(
    physical_group, tuple)):
        
        raise ValueError("The physical group as key of the constitutiv"+
        "e models dictionary must be either an integer or a tuple (for"+
        " multiple physical groups with the same constitutive model).")
    
    # Verifies if this or these physical groups are valid physical groups
    
    elif isinstance(physical_group, tuple):

        # Initializes a new physical group tuple

        sub_physicalGroupsList = []

        # Iterates through the physical groups in the tuple

        for group in physical_group:

            # If the key of the dictionary is a string, it is the physi-
            # cal group's name. Hence, converts it to its corresponding 
            # number tag

            group_tag = copy.deepcopy(group)

            print("Transforms the physical group name '"+str(group_tag)+
            "' to the tag:")
            
            if isinstance(group, str):

                try:

                    group_tag = physical_groupsNamesToTags[group]

                except:

                    raise KeyError("The physical group name '"+group+
                    "' was used in a variational form, but it does not"+
                    " exist in the dictionary of physical groups' name"+
                    "s to tags. This dictionary has the following keys"+
                    " and values: "+str(physical_groupsNamesToTags))
            
            print(group_tag, "\n")

            if not (group_tag in physical_groupsList):

                raise NameError("The physical group tag "+str(group_tag)
                +" was used to build the hyperelastic internal work, b"+
                "ut it is not a valid physical group. Here is the list"+
                " of the valid physical groups:\n"+str(
                physical_groupsList))
            
            sub_physicalGroupsList.append(group_tag)

        # Converts the physical group to this list

        physical_group = copy.deepcopy(sub_physicalGroupsList)
            
    else:

        if (not (physical_group in physical_groupsList)) and throw_error:

            raise NameError("The physical group tag "+str(physical_group
            )+" was used to build the hyperelastic internal work, but "+
            "it is not a valid physical group. Here is the list of the"+
            " valid physical groups:\n"+str(physical_groupsList))
        
    return physical_group

# Defines a function to project a field that is piecewise continuous in-
# to a finite element space. The field is given by a list of lists, whe-
# re each sublist is a pair of the function to be projected (i.e. the 
# value of the field) and the tag associated with its particular domain
# region (i.e. the physical group number or name)

@programming_tools.optional_argumentsInitializer({'solution_names': 
lambda: ["field", "field"]})

def project_piecewiseField(field_list, dx, V, physical_groupsList,
physical_groupsNamesToTags, solution_names=None, verbose=True):
    
    # Creates the projected field

    projected_fieldTrial = TrialFunction(V)

    v = TestFunction(V)

    projected_field = Function(V)

    projected_field.rename(*solution_names)

    # Assembles and solves the variational form

    bilinear_form = inner(projected_fieldTrial, v)*dx

    # Initializes the linear form

    linear_form = 0.0

    # Iterates through the pairs of function to be projected and physi-
    # cal group tag

    for projected_function, subdomain in field_list:

        if verbose:

            print("Projects over the subdomain "+str(subdomain))

        # Checks the subdomain for strings

        subdomain = verify_physicalGroups(subdomain,physical_groupsList, 
        physical_groupsNamesToTags=physical_groupsNamesToTags)

        # Updates the linear form

        linear_form += (inner(projected_function, v)*dx(subdomain))

    # Solves the algebraic system to get the FEM parameters

    solve(bilinear_form==linear_form, projected_field)

    # Returns the projected field

    return projected_field

# Defines a function to project a field over a region of the domain

def projection_overRegion(field, V, dx, subdomain, physical_groupsList,
physical_groupsNamesToTags):

    # Checks the subdomain for strings

    subdomain = verify_physicalGroups(subdomain, physical_groupsList, 
    physical_groupsNamesToTags=physical_groupsNamesToTags)

    # Creates the projected field

    projected_fieldTrial = TrialFunction(V)

    v = TestFunction(V)

    projected_field = Function(V)

    # Assembles and solves the variational form

    bilinear_form = inner(projected_fieldTrial, v)*dx
                     
    linear_form = (inner(field, v)*dx(subdomain))

    solve(bilinear_form==linear_form, projected_field)

    # Returns the projected field

    return projected_field