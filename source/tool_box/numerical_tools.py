# Routine to store functions for numerical analysis and numerical work-
# arounds

from dolfin import *

import ufl_legacy as ufl

########################################################################
#  Safe numerical operations to avoid lack of differentiability or di- #
#                            vision by zero                            #
########################################################################

# Defines a "safe" square root function, in order to avoid division by 
# zero in the variational form

def safe_sqrt(a):

    return ufl.sqrt(a+1.0e-15)

# Defines a function to give the closest value of a list to a given 
# float number

def closest_numberInList(number, list_ofNumbers):

    return min(list_ofNumbers, key=lambda x: abs(x-number))