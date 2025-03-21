# Routine to store some methods to help evaluate constitutive models. It
# includes stress tensors transformations and derivation. Cauchy and Mi-
# cropolar continua are considered

from dolfin import *

import ufl_legacy as ufl

import source.tool_box.tensor_tools as tensor_tools

# Defines the indices for Einstein summation notation

i, j, k, l = ufl.indices(4)

########################################################################
#                    Cauchy-continuum stress tensors                   #
########################################################################

# Defines a function to evaluate the second Piola-Kirchhoff stress ten-
# sor from the derivative of the Helmholtz potential w.r.t. the right
# Cauchy-Green strain tensor

def S_fromDPsiDC(helmholtz_potential, u):

    # Evaluates the deformation gradient

    I = Identity(3)

    F = grad(u)+I

    # Evaluates the right Cauchy-Green strain tensor, C. Makes C a 
    # variable to differentiate the Helmholtz potential with respect 
    # to C

    C = (F.T)*F
    
    C = variable(C)  

    # Evaluates the Helmholtz potential

    W = helmholtz_potential(C)

    # Evaluates the second Piola-Kirchhoff stress tensor differenti-
    # ating the potential w.r.t. C

    S = 2*diff(W,C)

    return S

# Defines a function to evaluate the Cauchy stress tensor as the push 
# forward of the second Piola-Kirchhoff stress tensor

def push_forwardS(S, u):

    # Evaluates the deformation gradient

    I = Identity(3)

    F = grad(u)+I

    # Evaluates the determinant of the deformation gradient

    J = ufl.det(F)

    # Pushes forward to the deformed configuration

    sigma = (1.0/J)*F*S*F.T

    return sigma

# Defines a function to transform the Cauchy stress to the second Piola-
# Kirchhoff stress

def S_fromCauchy(sigma, u):

    # Evaluates the deformation gradient

    I = Identity(3)

    F = grad(u)+I

    # Evaluates the determinant of the deformation gradient

    J = ufl.det(F)

    # Uses the pull back operation

    S = J*inv(F)*sigma*(inv(F).T)

    return S

# Defines a function to transform the Cauchy stress to the first Piola-
# Kirchhoff stress

def P_fromCauchy(sigma, u):

    # Evaluates the deformation gradient

    I = Identity(3)

    F = grad(u)+I

    # Evaluates the determinant of the deformation gradient

    J = ufl.det(F)

    # Uses the Piola transformation

    P = J*sigma*(inv(F).T)

    return P

# Defines a function to transform the first Piola-Kirchhoff stress from
# the second one

def P_fromS(S, u):

    # Evaluates the deformation gradient

    I = Identity(3)

    F = grad(u)+I

    # Transforms S to P

    P = F*S

    return P

########################################################################
#                     Micropolar curvature tensor                      #
########################################################################

# Defines the micropolar curvature tensor

def micropolar_curvatureTensor(phi):

    # Defines the permutation tensor for 3D euclidean space

    perm = ufl.PermutationSymbol(3)
    
    # Computes the micro-rotation tensor Rbar
     
    R_bar = tensor_tools.rotation_tensorEulerRodrigues(phi)
    
    # Computes the gradient of each component of R_bar to get a third 
    # order tensor, grad_Rbar

    grad_Rbar = grad(R_bar)

    # Constructs the third order curvature tensor K_third

    K_third = as_tensor(R_bar[j,i]*grad_Rbar[j,k,l], (i,k,l))

    # Contracts with the permutation tensor to form K_second

    K_second = 0.5*as_tensor(perm[i,j,k]*K_third[k,j,l], (i,l))

    return K_second