# Routine to store variational formulation methods

from dolfin import *

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
constitutive_modelDictionary, dx):
    
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

            verify_physicalGroups(physical_group, physical_groupsList)

            # Initializes objects for the stresses at the reference 
            # configuration

            first_piola = constitutive_model.first_piolaStress(
            trial_function)

            # Constructs the variational forms for the inner work

            inner_work += (inner(first_piola, grad(test_function))*
            dx(physical_group))

    # If the constitutive model is not a dictionary, the domain is homo-
    # geneous

    else:

        # Initializes objects for the stresses at the reference configu-
        # ration

        first_piola = constitutive_modelDictionary.first_piolaStress(
        trial_function)

        # Constructs the variational forms for the inner work

        inner_work = inner(first_piola, grad(test_function))*dx

    # Returns the inner work variational form

    return inner_work

########################################################################
#               Work done by Neumann boundary conditions               #
########################################################################

# Defines a function to construct the variational form of the work done
# by the traction vector in the reference configuration given a dictio-
# nary of traction loads, where the keys are the corresponding boundary
# physical groups and the values are the traction loads

def traction_work(traction_dictionary, field_variation, ds):

    # Gets the physical groups tags

    physical_groupsTags = set(ds.subdomain_data().array())

    # Initializes the variational form

    traction_form = 0.0

    # Iterates through the dictionary

    for physical_group, traction in traction_dictionary.items():

        # Verifies if this physical group is indeed in ds

        if not (physical_group in physical_groupsTags):

            raise NameError("The boundary physical group with tag "+str(
            physical_group)+" was attempted to construct the variation"+
            "al form of the traction work, but it does not exist insid"+
            "e the ds object. Probably the mesh does not have this bou"+
            "ndary physical group")

        traction_form += dot(traction, field_variation)*ds(
        physical_group)

    # Returns the variational form

    return traction_form

########################################################################
#                              Utilities                               #
########################################################################

# Defines a function to verify if a physical group is consistent and if
# it is the whole list of existing physical groups

def verify_physicalGroups(physical_group, physical_groupsList):

    # Tuples can be used as physical groups to integrate over multiple 
    # physical groups simultaneously

    if (not isinstance(physical_group, int)) and (not isinstance(
    physical_group, tuple)):
        
        raise ValueError("The physical group as key of the constitutiv"+
        "e models dictionary must be either an integer or a tuple (for"+
        " multiple physical groups with the same constitutive model).")
    
    # Verifies if this or these physical groups are valid physical groups
    
    elif isinstance(physical_group, tuple):

        # Iterates through the physical groups in the tuple

        for group in physical_group:

            if not (group in physical_groupsList):

                raise NameError("The physical group tag "+str(group)+
                " was used to build the hyperelastic internal work, bu"+
                "t it is not a valid physical group. Here is the list "+
                "of the valid physical groups:\n"+str(
                physical_groupsList))
            
    else:

        if not (physical_group in physical_groupsList):

            raise NameError("The physical group tag "+str(physical_group
            )+" was used to build the hyperelastic internal work, but "+
            "it is not a valid physical group. Here is the list of the"+
            " valid physical groups:\n"+str(physical_groupsList))

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