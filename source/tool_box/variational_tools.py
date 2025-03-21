# Routine to store variational formulation methods

from dolfin import *

# Defines a function to construct the variational form of a non-dissipa-
# tive solid mechanics problem. The constitutive model can be either a
# class (when the whole domain has only one constitutive model); or it 
# can be a dictionary, when the domain is heterogeneous. The keys of the
# dictionaries are the volumetric physical groups, whereas the values 
# are the constitutive model classes

def hyperelastic_internalWork(trial_function, test_function, 
constitutive_modelDictionary, dx):
    
    # Initializes the second order identity tensor and the deformation 
    # gradient
    
    I = Identity(3)

    F = grad(trial_function)+I

    # Initializes the variational form of the inner work

    inner_work = 0.0

    # If the constitutive model is a dictionary, the domain is heteroge-
    # neous

    if isinstance(constitutive_modelDictionary, dict):

        # Iterates through the dictionary

        for physical_group, constitutive_model in (
        constitutive_modelDictionary.items()):
            
            # Tuples can be used as physical groups to integrate over 
            # mutliple physical groups simultaneously

            if (not isinstance(physical_group, int)) and (not isinstance(
            physical_group, tuple)):
                
                raise ValueError("The physical group as key of the con"+
                "stitutive models dictionary must be either an integer"+
                " or a tuple (for multiple physical groups with the sa"+
                "me constitutive model).")

            # Initializes objects for the stresses at the reference 
            # configuration

            first_piola = constitutive_model.first_piolaStress(F)

            # Constructs the variational forms for the inner work

            inner_work += (inner(first_piola, grad(test_function))*
            dx(physical_group))

    # If the constitutive model is not a dictionary, the domain is homo-
    # geneous

    else:

        # Initializes objects for the stresses at the reference configu-
        # ration

        first_piola = constitutive_modelDictionary.first_piolaStress(F)

        # Constructs the variational forms for the inner work

        inner_work = inner(first_piola, grad(test_function))*dx

    # Returns the inner work variational form

    return inner_work

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