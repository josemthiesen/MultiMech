# Routine to store functions for numerical analysis and numerical work-
# arounds

from dolfin import *

import ufl_legacy as ufl

########################################################################
#                     Parametric loading functions                     #
########################################################################

# Defines a function to give different options of parametric loading 
# curves using ufl functions (so that they can be used in the variatio-
# nal ecosystem). These loading functions must spit out numbers between
# 0 and 1

def generate_loadingParametricCurves(curve_name):

    # Tests if it is linear

    if curve_name=="linear":

        return lambda x: x

    # Tests if it is the square root

    elif curve_name=="square root":

        return lambda x: ufl.sqrt(x)

    else:

        raise NameError("The parametric load curve '"+str(curve_name)+
        "' has not yet been implemented")

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