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

import source.tool_box.constitutive_tools as constitutive_tools

# Defines an abstract class to force all classes ahead to have the same
# methods. To enforce it, the abstract method is used before the methods.
# The material model classes have four methods in the case of hyperelas-
# ticity: the Helmholtz potential (strain energy); the second Piola-
# Kirchhoff stress tensor; the first Piola-Kirchhoff stress tensor; and
# the Cauchy stress tensor.

class HyperelasticMaterialModel(ABC):

    @abstractmethod

    # The following methods have the pass argument only because they 
    # will be defined in the child classes

    def strain_energy(self, strain_tensor):

        pass

    def second_piolaStress(self, displacement):

        pass

    def first_piolaStress(self, displacement):

        pass

    def cauchy_stress(self, displacement):

        pass

# Defines a class to evaluate the strain energy and the respective
# first Piola-Kirchhoff stress tensor for the neo-hookean hyperelastic 
# model from Javier Bonet, the same as in https://help.febio.org/docs/FE
# BioUser-3-6/UM36-4.1.3.17.html. It is called there as unconstrained
# Neo-Hookean material

class Neo_Hookean(HyperelasticMaterialModel):

    def __init__(self, material_properties):

        self.E = Constant(material_properties["E"])

        self.v = Constant(material_properties["v"])

        # Evaluates the Lamé parameters

        self.mu = Constant(self.E/(2*(1+self.v)))

        self.lmbda = Constant(self.v*self.E/((1+self.v)*(1-2*self.v)))

    # Defines a function to evaluate the Helmholtz free energy density

    def strain_energy(self, C):

        # Evaluates the right Cauchy-Green strain tensor invariants

        I1_C = ufl.tr(C)

        I2_C = ufl.det(C)

        J = ufl.sqrt(I2_C)
        
        # Evaluates the trace-related part

        psi_1 = (self.mu/2)*(I1_C - 3)

        # Evaluates the jacobian-related part

        psi_2 = -(self.mu*ufl.ln(J))+((self.lmbda*0.5)*((ufl.ln(J))**2))

        return psi_1+psi_2
    
    # Defines a function to evaluate the second Piola-Kirchhoff stress 
    # tensor as the derivative of the Helmholtz free energy density po-
    # tential

    def second_piolaStress(self, u):

        S = constitutive_tools.S_fromDPsiDC(self.strain_energy, u)

        return S

    # Defines a function to evaluate the first Piola-Kirchhoff stress 
    # tensor as the result of the proper operation over the second one

    def first_piolaStress(self, u):

        # Evaluates the second Piola-Kirchhoff stress tensor and pulls
        # it back to the first Piola-Kirchhoff stress tensor

        S = self.second_piolaStress(u)
        
        P = constitutive_tools.P_fromS(S, u)

        return P
    
    # Defines a function to evaluate the Cauchy stress tensor as the 
    # push forward of the second Piola-Kirchhoff stress tensor
    
    def cauchy_stress(self, u):

        # Evaluates the second Piola-Kirchhoff stress tensor and pushes 
        # it forward to the deformed configuration

        S = self.second_piolaStress(u)

        sigma = constitutive_tools.push_forwardS(S, u)

        return sigma