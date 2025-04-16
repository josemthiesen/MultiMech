# Routine to store variational formulation methods

from dolfin import *

import copy

import source.tool_box.tensor_tools as tensor_tools

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
constitutive_modelDictionary, dx, domain_physGroupsNamesToTags=dict(),
verbose=False):
    
    # Gets the physical groups from the domain mesh function

    physical_groupsList = set(dx.subdomain_data().array())

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
            domain_physGroupsNamesToTags)

            # Initializes objects for the stresses at the reference 
            # configuration

            first_piola = constitutive_model.first_piolaStress(
            trial_function)

            # If the physical group is indeed a list of tags, i.e. mul-
            # tiple regions are integrated using the same constitutive
            # model

            if isinstance(physical_group, list):

                # Iterates through the subdomains

                for sub_physicalGroup in physical_group:

                    print("The volume of the", sub_physicalGroup, "sub"+
                    "domain is:", assemble(1*dx(sub_physicalGroup)), 
                    "\n")

                    # Constructs the variational forms for the inner 
                    # work

                    inner_work += (inner(first_piola, grad(test_function
                    ))*dx(sub_physicalGroup))

            else:

                print("The volume of the", physical_group, "subdomain "+
                "is:", assemble(1*dx(physical_group)), "\n")

                # Constructs the variational forms for the inner work

                inner_work += (inner(first_piola, grad(test_function))*
                dx(physical_group))

    # If the constitutive model is not a dictionary, the domain is homo-
    # geneous

    else:

        print("The volume of the domain is:", assemble(1*dx), "\n")

        # Initializes objects for the stresses at the reference configu-
        # ration

        first_piola = constitutive_modelDictionary.first_piolaStress(
        trial_function)

        # Constructs the variational forms for the inner work

        inner_work = inner(first_piola, grad(test_function))*dx

    # Returns the inner work variational form

    if verbose:

        print("Finishes creating the variational form of the inner wor"+
        "k done by the first Piola stress tensor in a Cauchy continuum"+
        " medium\n")

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
constitutive_modelDictionary, dx, domain_physGroupsNamesToTags=dict(),
verbose=False):
    
    # Gets the physical groups from the domain mesh function

    physical_groupsList = set(dx.subdomain_data().array())

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
            domain_physGroupsNamesToTags)

            # Initializes objects for the stresses at the reference 
            # configuration

            first_piola, couple_firstPiola = constitutive_model.first_piolaStress(
            [displacement_trialFunction, microrotation_trialFunction])

            kirchhoff, couple_kirchhoff = constitutive_model.kirchhoff_stress(
            [displacement_trialFunction, microrotation_trialFunction])

            # If the physical group is indeed a list of tags, i.e. mul-
            # tiple regions are integrated using the same constitutive
            # model

            if isinstance(physical_group, list):

                # Iterates through the individual physical groups

                for sub_physicalGroup in physical_group:

                    # Constructs the variational forms for the inner 
                    # work of the first Piola-Kirchhoff stress

                    inner_work += (inner(first_piola, grad(
                    displacement_testFunction))*dx(sub_physicalGroup))

                    # Adds the parcel of the couple stress

                    inner_work += ((inner(couple_firstPiola, grad(
                    microrotation_testFunction))*dx(sub_physicalGroup))-
                    (inner(kirchhoff, tensor_tools.skew_2OrderTensor(
                    microrotation_testFunction))*dx(sub_physicalGroup)))

            else:

                # Constructs the variational forms for the inner work of
                # the first Piola-Kirchhoff stress

                inner_work += (inner(first_piola, grad(
                displacement_testFunction))*dx(physical_group))

                # Adds the parcel of the couple stress

                inner_work += ((inner(couple_firstPiola, grad(
                microrotation_testFunction))*dx(physical_group))-(inner(
                kirchhoff, tensor_tools.skew_2OrderTensor(
                microrotation_testFunction))*dx(physical_group)))

    # If the constitutive model is not a dictionary, the domain is homo-
    # geneous

    else:

        # Initializes objects for the stresses at the reference configu-
        # ration

        first_piola, couple_firstPiola = constitutive_modelDictionary.first_piolaStress(
        [displacement_trialFunction, microrotation_trialFunction])

        kirchhoff, couple_kirchhoff = constitutive_modelDictionary.kirchhoff_stress(
        [displacement_trialFunction, microrotation_trialFunction])

        # Constructs the variational forms for the inner work of the
        # first Piola-Kirchhoff stress

        inner_work += (inner(first_piola, grad(displacement_testFunction
        ))*dx)

        # Adds the parcel of the couple stress

        inner_work += ((inner(couple_firstPiola, grad(
        microrotation_testFunction))*dx)-(inner(kirchhoff, 
        tensor_tools.skew_2OrderTensor(microrotation_testFunction))*dx))

    # Returns the inner work variational form

    if verbose:

        print("Finishes creating the variational form of the inner wor"+
        "k done by the first Piola stress tensor and its couple stress"+
        " in a micropolar continuum medium\n")

    return inner_work

########################################################################
#               Work done by Neumann boundary conditions               #
########################################################################

# Defines a function to construct the variational form of the work done
# by the traction vector in the reference configuration given a dictio-
# nary of traction loads, where the keys are the corresponding boundary
# physical groups and the values are the traction loads

def traction_work(traction_dictionary, field_variation, ds, 
boundary_physGroupsNamesToTags=dict(), verbose=False):

    # Gets the physical groups tags

    physical_groupsTags = set(ds.subdomain_data().array())

    # Initializes the variational form

    traction_form = 0.0

    # Iterates through the dictionary

    for physical_group, traction in traction_dictionary.items():

        # Verifies if this physical group is indeed in ds
            
        physical_group = verify_physicalGroups(physical_group, 
        physical_groupsTags, physical_groupsNamesToTags=
        boundary_physGroupsNamesToTags)

        if isinstance(physical_group, list):

            for sub_physicalGroup in physical_group:

                print("The physical group "+str(sub_physicalGroup)+" h"+
                "as an area of "+str(assemble(1*ds(sub_physicalGroup)))+
                "\n")

                traction_form += dot(traction, field_variation)*ds(
                sub_physicalGroup)

        else:

            print("The physical group "+str(physical_group)+" has an a"+
            "rea of "+str(assemble(1*ds(physical_group)))+"\n")

            traction_form += dot(traction, field_variation)*ds(
            physical_group)

    # Returns the variational form

    if verbose:

        print("Finishes creating the variational form of the work done"+
        " by the traction on the boundary\n")

    return traction_form

########################################################################
#                              Utilities                               #
########################################################################

# Defines a function to verify if a physical group is consistent and if
# it is the whole list of existing physical groups

def verify_physicalGroups(physical_group, physical_groupsList, 
physical_groupsNamesToTags=dict()):
            
    # If the key of the dictionary is a string, it is the physical
    # group's name. Hence, converts it to its corresponding number tag
    
    if isinstance(physical_group, str):

        try:

            print("Transforms the physical group name '"+physical_group+
            "' to the tag:")

            physical_group = physical_groupsNamesToTags[physical_group]

            print(physical_group, "\n")

        except:

            raise KeyError("The physical group name '"+physical_group+
            "' was used in a variational form, but it does not exist i"+
            "n the dictionary of physical groups' names to tags. This "+
            "dictionary has the following keys and values: "+str(
            physical_groupsNamesToTags))

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

        if not (physical_group in physical_groupsList):

            raise NameError("The physical group tag "+str(physical_group
            )+" was used to build the hyperelastic internal work, but "+
            "it is not a valid physical group. Here is the list of the"+
            " valid physical groups:\n"+str(physical_groupsList))
        
    return physical_group

# Defines a function to project a field over a region of the domain

def projection_overRegion(field, V, dx, subdomain):

    # Creates the projected field

    projected_fieldTrial = TrialFunction(V)

    v = TestFunction(V)

    projected_field = Function(V)

    # Assembles and solves the variational form

    solve((inner(projected_fieldTrial, v)*dx==(inner(field, v)*dx(
    subdomain))), projected_field)

    # Returns the projected field

    return projected_field