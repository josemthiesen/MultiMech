# Routine to store variational formulation methods

from dolfin import *

# Defines a function to construct the variational form of the work done
# by the traction vector in the reference configuration given a dictio-
# nary of traction loads, where the keys are the corresponding boundary
# physical groups and the values are the traction loads

def traction_work(traction_dictionary, field_variation, ds):

    # Initializes the variational form

    traction_form = 0.0

    # Iterates through the dictionary

    for physical_group, traction in traction_dictionary.items():

        traction_form += dot(traction, field_variation)*ds(
        physical_group)

    # Returns the variational form

    return traction_form