# Routine to store tensorial operations of tensor algebra and tensor 
# calculus

from dolfin import *

import ufl_legacy as ufl

import numpy as np

import source.tool_box.numerical_tools as n_tools

# Define the indices for Einstein summation notation

i, j, k, l = ufl.indices(4)

########################################################################
#                                Norms                                 #
########################################################################

# Defines a function to evaluate the L2 norm of a vector

def L2_normVector(vector):

    return n_tools.safe_sqrt(dot(vector, vector))

########################################################################
#                        Skew tensor operations                        #
########################################################################

# Defines a function to construct the skew-symmetric second order tensor
# W from its axial vector

def skew_2OrderTensor(axial_vector):
   
   # Evaluates the permutation third order tensor
   
   perm = ufl.PermutationSymbol(3)

   # Constructs the skew second order tensor W
   
   W = as_tensor(perm[j,i,k]*axial_vector[k], (i,j))

   return W

########################################################################
#                               Rotation                               #
########################################################################

# Defines a function to evaluate a rotation second order tensor using 
# the Euler-Rodrigues formula, given the corresponding axial vector, 
# which sets the axis where the rotation is carried about, as well as 
# the angle of rotation in radians from its L2 norm

def rotation_tensorEulerRodrigues(phi):

    # Evaluates the rotation angle

    rotation_angle = L2_normVector(phi)

    # Evaluates the skew tensor and the identity tensor

    I = Identity(3)

    W = skew_2OrderTensor(phi)

    # Evaluates the Euler-Rodrigues formula

    R_bar = ((cos(rotation_angle)*I)+((sin(rotation_angle)/
    rotation_angle)*W)+(((1-cos(rotation_angle))/(rotation_angle**2))*
    tensor_productVectorVector(phi,phi)))

    return R_bar

########################################################################
#                           Tensor products                            #
########################################################################

# Defines the tensor product between two vectors, a and b

def tensor_productVectorVector(a,b):

    C = ufl.as_tensor(a[i]*b[j], (i,j))

    return C