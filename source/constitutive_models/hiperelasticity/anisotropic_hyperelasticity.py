# Routine to store anisotropic constitutive models

from abc import ABC, abstractmethod

import numpy as np

from dolfin import *

import ufl_legacy as ufl

import source.tool_box.tensor_tools as tensor_tools

# Defines an abstract class to force all classes ahead to have the same
# methods. To enforce it, the abstract method is used before the methods

class MaterialModel(ABC):

    @abstractmethod

    # The following methods have the pass argument only because they 
    # will be defined in the child classes

    def strain_energy(self, Measure):

        pass

    def stress_tensor(self, Measure):

        pass

# Defines a function to implement the strain energy and the stress ten-
# sor for the Holzapfel-Gasser-Ogden unconstrained material model. This
# implementation was done as in https://help.febio.org/docs/FEBioUser-3-
# 6/UM36-4.1.4.11.html, where it is called unconstrained Holzapfel-
# Gasser-Ogden

class Holzapfel_Gasser_Ogden_Unconstrained(MaterialModel):

    # Initializes the properties

    def __init__(self, material_properties):
        
        self.c = material_properties["c"]

        self.k1 = material_properties["k1"]

        self.k2 = material_properties["k2"]

        self.gamma = material_properties["gamma"]*(ufl.pi/180.0)

        self.kappa = material_properties["kappa"]

        self.k = material_properties["k"]

        self.e1 = material_properties["local system of coordinates: a "+
        "direction"]

        self.e2 = material_properties["local system of coordinates: d "+
        "direction"]

        self.e3 = ufl.cross(self.e1, self.e2)

        # Normalizes these vectors

        norm_e1 = tensor_tools.L2_normVector(self.e1)

        norm_e2 = tensor_tools.L2_normVector(self.e2)

        norm_e3 = tensor_tools.L2_normVector(self.e3)

        self.e1 = (1/norm_e1)*self.e1

        self.e2 = (1/norm_e2)*self.e2

        self.e3 = (1/norm_e3)*self.e3

    # Defines the strain energy

    def strain_energy(self, C):

        # Evaluates the invariants of the right Cauhcy-Green tensor
        
        J = ufl.variable(ufl.sqrt(ufl.det(C)))

        I1_C = ufl.variable(ufl.tr(C))

        # Defines the rotation axis vector
        n = self.e3
        norm_n = tensor_tools.L2_normVector(self.e3)

        n = n/norm_n

        W = ufl.as_tensor([[0, -n[2], n[1]],
                       [n[2], 0, -n[0]],
                       [-n[1], n[0], 0]])

        W = ufl.variable(W)

        # Identity matrix
        I = ufl.Identity(3)
        I = ufl.variable(I)

        R1 = ufl.cos(self.gamma) * I + ufl.sin(self.gamma) * W + (1 - ufl.cos(self.gamma)) * (ufl.outer(n, n))
        R2 = ufl.cos(-self.gamma) * I + ufl.sin(-self.gamma) * W + (1 - ufl.cos(-self.gamma)) * (ufl.outer(n, n))

        alpha_1 = R1 * self.e1
        alpha_2 = R2 * self.e1

        I4_alpha_1 = ufl.dot(alpha_1, (C * alpha_1))
        I4_alpha_2 = ufl.dot(alpha_2, (C * alpha_2))

        I4_alpha_1 = ufl.variable(I4_alpha_1)
        I4_alpha_2 = ufl.variable(I4_alpha_2)

        def Macaulay(variable):
            return ufl.conditional(ufl.lt(variable, 0), 0, variable)

        E_alpha_1 = self.kappa * (I1_C - 3) + (1 - 3 * self.kappa) * (I4_alpha_1 - 1)
        E_alpha_1 = ufl.variable(E_alpha_1)
        E_alpha_2 = self.kappa * (I1_C - 3) + (1 - 3 * self.kappa) * (I4_alpha_2 - 1)
        E_alpha_2 = ufl.variable(E_alpha_2)

        energy_matrix = self.c * (I1_C - 3) 

        energy_fiber_1 = (self.k1 / (2 * self.k2)) * (ufl.exp(self.k2 * ((Macaulay(E_alpha_1) ** 2))) - 1)
        energy_fiber_2 = (self.k1 / (2 * self.k2)) * (ufl.exp(self.k2 * ((Macaulay(E_alpha_2) ** 2))) - 1)

        energy_volumetric = ((self.k) * ((((J**2)-1)/2) - ufl.ln(J)))  - (2*self.c * ufl.ln(J))

        return (energy_matrix + energy_fiber_1 + energy_fiber_2 + energy_volumetric)

    def stress_tensor(self, F):
        
        C = (F.T)*F
        
        C = ufl.variable(C)

        W = self.strain_energy(C)

        # Second Piola-Kirchhoff stress tensor
        S = 2 * diff(W, C)

        # First Piola-Kirchhoff stress tensor
        P = F * S

        return P

# Implementation of the neo-hookean hyperelastic model from Javier Bonet, the same as in 
# https://help.febio.org/docs/FEBioUser-3-6/UM36-4.1.3.17.html called unconstrained 
# Neo-Hookean material

class Neo_Hookean_Hyperelasticity(MaterialModel):

    def __init__(self, material_properties, localcsys_properties):
        self.E = material_properties["E"]
        self.v = material_properties["v"]
        self.mu = self.E/(2*(1+self.v))
        self.lambda_ = self.v*self.E/((1+self.v)*(1-2*self.v))

    def strain_energy(self, C):

        I1_C = tr(C)
        III_C = det(C)

        J = sqrt(III_C)
        
        # This is a compressible neo-Hookean material. It is derived from the following hyperelastic strain energy function
        energy_1 = (self.mu/2)*(I1_C - 3)
        energy_2 = -self.mu*ln(J)
        energy_3 = (self.lambda_/2)*(ln(J))**2

        return energy_1 + energy_2 + energy_3

    def stress_tensor(self, F):

        F = variable(F)
        C = (F.T)*F
        C = variable(C)

        J = sqrt(det(C))    
        W = self.strain_energy(C)

        # Second Piola-Kirchhoff stress tensor
        S = 2*diff(W,C)

        # First Piola-Kirchhoff stress tensor
        P = F*S

        return P