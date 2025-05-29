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

    # If the curve name is a list, gets parameters for the load curve

    parameters_curve = None

    if isinstance(curve_name, list):

        if len(curve_name)>1:

            parameters_curve = curve_name[1]

            curve_name = curve_name[0]

            # Tests if the parameters_curve is a dictionary

            if not isinstance(parameters_curve, dict):

                raise TypeError("You've set the parametric load curve "+
                "as a list, i.e. the first component is the name of th"+
                "e curve, and the second component is a dictionary of "+
                "optional parameters for the load curve. The problem i"+
                "s: you haven't set the optional parameters as a dicti"+
                "onary, rather as "+str(parameters_curve))

        elif len(curve_name)>0:

            curve_name = curve_name[0]

        else:

            raise IndexError("If a load parametric curve is given as a"+
            " list, it has to have two elements: the first one is the "+
            "curve's name, and the second one is the dictionary of par"+
            "ameters")

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