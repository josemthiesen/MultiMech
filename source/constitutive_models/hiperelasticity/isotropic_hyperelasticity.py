# Routine to store anisotropic constitutive models
# 
# The list of implemented isotropic hyperelastic models is:
# 
# 1. Neo-Hookean

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