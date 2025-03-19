from abc import ABC, abstractmethod
import numpy as np
from dolfin import *
import ufl_legacy as ufl

class MaterialModel(ABC):
    @abstractmethod
    def strain_energy(self, Measure):
        pass # This method attribute will be implemented by the subclasses
    def stress_tensor(self, Measure):
        pass # This method attribute will be implemented by the subclasses

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

        I_C = tr(C)
        III_C = det(C)

        J = sqrt(III_C)
        
        # This is a compressible neo-Hookean material. It is derived from the following hyperelastic strain energy function
        energy_1 = (self.mu/2)*(I_C - 3)
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

# Implementation of the Holzapfel-Gasser-Ogden Unconstrained model, the same as in 
# https://help.febio.org/docs/FEBioUser-3-6/UM36-4.1.4.11.html called unconstrained 
# Holzapfel-Gasser-Ogden

class Holzapfel_Gasser_Ogden_Unconstrained(MaterialModel):

    def __init__(self, material_properties, localcsys_properties):
        
        self.c = material_properties["c"]
        self.k1 = material_properties["k1"]
        self.k2 = material_properties["k2"]
        self.gamma = material_properties["gamma"]*ufl.pi/180.0
        self.kappa = material_properties["kappa"]
        self.k = material_properties["k"]
        self.e1 = localcsys_properties["a"]
        self.e2 = localcsys_properties["d"]
        self.e3 = ufl.cross(self.e1, self.e2)

    def strain_energy(self, C):

        J = ufl.sqrt(ufl.det(C))
        J = ufl.variable(J)

        I_C = ufl.tr(C)

        I_C = ufl.variable(I_C)

        # Define the unit rotation axis
        n = self.e3
        norm_n = sqrt(n[0]**2 + n[1]**2 + n[2]**2)
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

        E_alpha_1 = self.kappa * (I_C - 3) + (1 - 3 * self.kappa) * (I4_alpha_1 - 1)
        E_alpha_1 = ufl.variable(E_alpha_1)
        E_alpha_2 = self.kappa * (I_C - 3) + (1 - 3 * self.kappa) * (I4_alpha_2 - 1)
        E_alpha_2 = ufl.variable(E_alpha_2)

        energy_matrix = self.c * (I_C - 3) 

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