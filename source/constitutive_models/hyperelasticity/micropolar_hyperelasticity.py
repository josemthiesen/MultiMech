# Routine to store anisotropic constitutive models
# 
# The list of implemented isotropic hyperelastic models is:
# 
# 1. Neo-Hookean

from abc import ABC, abstractmethod

from dolfin import *

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

    def strain_energy(self, strain_tensor, curvature_tensor):

        pass

    def second_piolaStress(self, displacement, microrotation):

        pass

    def first_piolaStress(self, displacement, microrotation):

        pass

    def cauchy_stress(self, displacement, microrotation):

        pass

# Defines a class to evaluate the strain energy and the respective
# first Piola-Kirchhoff stress tensor for the neo-hookean hyperelastic 
# micropolar model. This implementation is the same as from Bauer et al-
# Micropolar hyperelasticity: constitutive model, consistent lineariza-
# tion and simulation of 3D scale effects

class Micropolar_Neo_Hookean(HyperelasticMaterialModel):

    def __init__(self, material_properties):

        self.mu = material_properties["mu"]

        self.lmbda = material_properties["lambda"]

        self.kappa = material_properties["kappa"]

        self.alpha = material_properties["alpha"]

        self.beta = material_properties["beta"]

        self.gamma = material_properties["gamma"]

    # Defines a function to evaluate the Helmholtz free energy density

    def strain_energy(self, V_bar, k_curvatureSpatial):

        # Evaluates the traces

        I1 = tr(V_bar*V_bar.T)

        I2 = tr(V_bar*V_bar)

        I3 = tr(k_curvatureSpatial)

        I4 = tr(k_curvatureSpatial*k_curvatureSpatial)

        I5 = tr(k_curvatureSpatial*k_curvatureSpatial.T)

        # And the determinant

        J = det(V_bar)

        # Evaluates the energy parcels from Bauer et al

        psi_NH = 0.5*self.mu*(I1-3.0)

        psi_vol = ((0.25*self.lmbda*((J**2)-1))-(0.5*self.lmbda*(ln(J)))
        -self.mu*ln(J))

        psi_hat = 0.25*self.kappa*(I1-I2)

        psi_tilde = (0.5*((self.alpha*(I3**2))+(self.beta*I4)+(
        self.gamma*I5)))

        return psi_NH+psi_vol+psi_hat+psi_tilde
    
    # Defines a function to evaluate the second Piola-Kirchhoff stress 
    # tensors using a pull-back operation over the Cauchy stress tensors.
    # The argument is a list of the two fields

    def second_piolaStress(self, fields_list):

        # Retrieves the fields

        u, phi = fields_list

        # Evaluates the Cauchy stress and the couple stress

        sigma, sigma_couple = self.cauchy_stress(fields_list)

        # Evaluates the second Piola-Kirchhoff stress tensors using the
        # pull back operation

        S = constitutive_tools.S_fromCauchy(sigma, u)

        S_couple = constitutive_tools.S_fromCauchy(sigma_couple, u)

        return S, S_couple

    # Defines a function to evaluate the first Piola-Kirchhoff stress 
    # and its couple from the Cauchy stress

    def first_piolaStress(self, fields_list):

        # Retrieves the fields

        u, phi = fields_list

        # Evaluates the Cauchy stress and the couple stress

        sigma, sigma_couple = self.cauchy_stress(fields_list)

        # Evaluates the first Piola-Kirchhoff stress tensors using the
        # Piola transformation

        P = constitutive_tools.P_fromCauchy(sigma, u)

        P_couple = constitutive_tools.P_fromCauchy(sigma_couple, u)

        return P, P_couple

    # Defines a function to evaluate the Kirchhoff stress and its couple 
    # from the Cauchy stress

    def kirchhoff_stress(self, fields_list):

        # Retrieves the fields

        u, phi = fields_list

        # Evaluates the Cauchy stress and the couple stress

        sigma, sigma_couple = self.cauchy_stress(fields_list)

        # Multiplies by the jacobian to get the Kirchhoff stresses

        tau = constitutive_tools.tau_fromCauchy(sigma, u)

        tau_couple = constitutive_tools.tau_fromCauchy(sigma_couple, u)

        return tau, tau_couple
    
    # Defines a function to evaluate the Cauchy stress tensor using the 
    # derivative of the Helmholtz potential with respect to the V_bar 
    # tensor. This function receives the deformation gradient, the rota-
    # tion tensor of the microrotation field, and the curvature tensor 
    # in the referential configuration
    
    def cauchy_stress(self, fields_list):

        # Retrieves the fields

        u, phi = fields_list

        # Evaluates the deformation gradient

        I = Identity(3)

        F = grad(u)+I 

        # Evaluates the rotation tensor using phi

        R_bar = tensor_tools.rotation_tensorEulerRodrigues(phi)

        # Evaluates the curvature tensor in the referential configuration

        K_curvatureReferential = constitutive_tools.micropolar_curvatureTensor(
        phi)

        # Evaluates the micropolar stretch and the jacobian
    
        V_bar = F*(R_bar.T)

        # Evaluates the push-forward of the curvature tensor

        k_curvatureSpatial = R_bar*K_curvatureReferential*R_bar.T

        # Transforms the tensors into variables to differentiate the e-
        # nergy potential

        #V_barVar = variable(V_bar)

        #k_curvatureVar = variable(k_curvatureSpatial)

        # Evaluates the total energy density

        #psi_total1 = self.strain_energy(V_barVar, k_curvatureSpatial)

        #dPsi1 = diff(psi_total1, V_barVar)

        #sigma = V_bar*dPsi1.T

        #psi_total2 = self.strain_energy(V_bar, k_curvatureVar)

        #dPsi2 = diff(psi_total2, k_curvatureVar)

        #sigma_couple = V_bar*dPsi2.T

        # Evaluates the Cauchy and the couple Cauchy stress tensors

        """

        k_curvatureSpatialTransposed = variable(k_curvatureSpatial.T)

        V_barTransposed = variable(V_bar.T)

        psi_total = self.strain_energy(V_barTransposed.T, 
        k_curvatureSpatialTransposed.T)

        sigma = V_bar*diff(psi_total, V_barTransposed)

        sigma_couple = V_bar*diff(psi_total,k_curvatureSpatialTransposed)"""

        #"""
        J = det(V_bar)

        sigma = (((self.lmbda/2)*((J*J)-1)*I)+(self.mu*((V_bar*V_bar.T)-
        I))+((self.kappa/2)*((V_bar*V_bar.T)-(V_bar*V_bar))))

        sigma_couple = V_bar*((self.alpha*tr(k_curvatureSpatial)*I)+(
        self.beta*k_curvatureSpatial)+(self.gamma*k_curvatureSpatial.T))#"""

        # Returns them

        return sigma, sigma_couple