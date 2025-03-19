# Routine to store variational formulation methods

from dolfin import *

# Defines a function to construct the variational form of the work done
# by the traction vector in the reference configuration given a dictio-
# nary of traction loads, where the keys are the corresponding boundary
# physical groups and the values are the traction loads

def traction_work(traction_dictionary, field_variation, ds):

    # Gets the physical groups tags

    physical_groupsTags = set(ds.sub_domain_data().array())

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