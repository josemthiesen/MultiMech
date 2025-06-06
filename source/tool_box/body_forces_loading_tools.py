# Routine to store tools to create loads

from dolfin import *

import ufl_legacy as ufl

import numpy as np

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.tensor_tools as tensor_tools

import source.tool_box.numerical_tools as numerical_tools

########################################################################
#                          Referential forces                          #
########################################################################

# Defines a function to construct a uniform body force on a domain in 
# the referential configuration

def UniformReferentialTraction(amplitude_bodyX, amplitude_bodyY,
amplitude_bodyZ, t=0.0, t_final=1.0, parametric_load_curve="linear"):
    
    # Verifies if the parametric load curve is a string and, then, uses
    # it to convert to an actual function

    if isinstance(parametric_load_curve, str):

        parametric_load_curve = numerical_tools.generate_loadingParametricCurves(
        parametric_load_curve)

    # Creates the time constants

    time_constant = Constant(t)

    maximum_time = Constant(t_final)

    # Evaluates the body force vector at the referential configuration

    B = as_vector([parametric_load_curve(time_constant/maximum_time)*
    amplitude_bodyX, parametric_load_curve(time_constant/
    maximum_time)*amplitude_bodyY, parametric_load_curve(
    time_constant/maximum_time)*amplitude_bodyZ])

    # Returns the body force vector and the time constant

    return B, time_constant

########################################################################
#                           Follower forces                            #
########################################################################