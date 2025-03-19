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