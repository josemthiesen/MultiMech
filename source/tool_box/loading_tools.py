# Routine to store tools to create loads

from dolfin import *

########################################################################
#                           Follower forces                            #
########################################################################

# Defines a function to construct a normal uniform traction in a surfa-
# ce, this means the traction vector is always perpendicular to the sur-
# face. The user can set how the amplitude of the function is interpola-
# ted between the initial and final time points; the default value is a
# linear curve

def normal_uniformForce(mesh_dataClass, displacement_field, 
maximum_traction, t=0.0, t_max=1.0, amplitude_curve=lambda x: x):

    # Creates the time constants

    time_constant = Constant(t)

    maximum_time = Constant(t_max)

    # Evaluates the deformation gradient

    I = Identity(3)

    F = grad(displacement_field)+I

    J = det(F)

    # Evaluates the traction vector at the referential configuration

    T = (amplitude_curve(time_constant/maximum_time)*maximum_traction*J*
    ((inv(F)).T)*mesh_dataClass.n)

    # Returns the traction and the time constant

    return T, time_constant