# Routine to store tools to create loads

from dolfin import *

########################################################################
#                          Referential forces                          #
########################################################################

# Defines a function to construct a uniform traction on a surface in the
# referential configuration

def UniformReferentialTraction(maximum_tractionX, maximum_tractionY,
maximum_tractionZ, t=0.0, t_max=1.0, amplitude_curve=lambda x: x):

    # Creates the time constants

    time_constant = Constant(t)

    maximum_time = Constant(t_max)

    # Evaluates the traction vector at the referential configuration

    T = as_vector([amplitude_curve(time_constant/maximum_time)*
    maximum_tractionX, amplitude_curve(time_constant/maximum_time)*
    maximum_tractionY, amplitude_curve(time_constant/maximum_time)*
    maximum_tractionZ])

    print("Cria")

    # Returns the traction and the time constant

    return T, time_constant

########################################################################
#                           Follower forces                            #
########################################################################

# Defines a function to construct a normal uniform traction in a surfa-
# ce, this means the traction vector is always perpendicular to the sur-
# face. The user can set how the amplitude of the function is interpola-
# ted between the initial and final time points; the default value is a
# linear curve

def NormalUniformFollowerTraction(field, mesh_dataClass,
maximum_traction, t=0.0, t_max=1.0, amplitude_curve=lambda x: x):

    # Creates the time constants

    time_constant = Constant(t)

    maximum_time = Constant(t_max)

    # Evaluates the deformation gradient

    I = Identity(3)

    F = grad(field)+I

    J = det(F)

    # Evaluates the traction vector at the referential configuration

    T = (amplitude_curve(time_constant/maximum_time)*maximum_traction*J*
    ((inv(F)).T)*mesh_dataClass.n)

    # Returns the traction and the time constant

    return T, time_constant