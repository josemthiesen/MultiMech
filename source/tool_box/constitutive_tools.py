# Routine to store some methods to help evaluate constitutive models

from dolfin import *

import ufl_legacy as ufl

# Defines a function to evaluate the second Piola-Kirchhoff stress ten-
# sor from the derivative of the Helmholtz potential w.r.t. the right
# Cauchy-Green strain tensor

def S_fromDPsiDC(F, helmholtz_potential):

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

def push_forwardS(S, F):

    # Evaluates the determinant of the deformation gradient

    J = ufl.det(F)

    # Pushes forward to the deformed configuration

    sigma = (1.0/J)*F*S*F.T

    return sigma