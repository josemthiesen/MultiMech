from dolfin import *

from mpi4py import MPI

import ufl_legacy as ufl

import numpy as np

import matplotlib.pyplot as plt

#import periodic_structure as mesher

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.numerical_tools as n_tools

import source.tool_box.tensor_tools as t_tools

import source.tool_box.functional_tools as fun_tools

# Define the indices for Einstein summation notation

i, j, k, l = ufl.indices(4)

########################################################################
########################################################################
##                      User defined parameters                       ##
########################################################################
########################################################################

########################################################################
#                         Material properties                          #
########################################################################

# Sets a dictionary for material property for each physical group in the
# volumetric domain. The keys are the physical groups and the values are
# the material properties

mu_materials = dict()

mu_materials[3] = 26.12

mu_materials[4] = 26.12

# Defines some micropolar constitutive parameters

alpha = 0.0

beta = 0.0

gamma = 5.22e-8

########################################################################
#                          Macro deformations                          #
########################################################################

# A set of macro deformations are stored in txt files. One must define
# their file names

macro_displacementName = "u_RVE_homogenized"

macro_gradDisplacementName = "grad_u_RVE_homogenized"

macro_microrotationName = "phi_RVE_homogenized"

macro_gradMicrorotationName = "grad_phi_RVE_homogenized"

########################################################################
#                                 Mesh                                 #
########################################################################

# Defines the name of the file to save the mesh in. Do not write the fi-
# le termination, e.g. .msh or .xdmf; both options will be saved automa-
# tically

file_name = "test_meshes//RVE_mesh"

# Defines a flag to generate a new mesh or not

flag_newMesh = True

########################################################################
#                            Function space                            #
########################################################################

# Defines the shape functions degree

polynomial_degree = 1

########################################################################
#                           Solver parameters                          #
########################################################################

# Sets some parameters

parameters["form_compiler"]["representation"] = "uflacs"
parameters["allow_extrapolation"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 2

# Sets the solver parameters

linear_solver = "minres"

relative_tolerance = 1e-2

maximum_iterations = 50

preconditioner = "petsc_amg"

krylov_absoluteTolerance = 1e-6

krylov_relativeTolerance = 1e-6

krylov_maximumIterations = 15000

krylov_monitorConvergence = False

# Sets the final pseudotime of the simulation

t_final = 1.0

# Sets the maximum number of steps of loading

maximum_loadingSteps = 11

########################################################################
########################################################################
##               Calculation: Mr. User, take care ahead!              ##
########################################################################
########################################################################

########################################################################
#                                 Mesh                                 #
########################################################################

# Generates the mesh and writes it 

#if flag_newMesh:

#    mesher.generate_periodicMesh(file_name, flag_transfinite=0, verbose=
#    True, lc=5e-4)

# Reads the mesh and constructs some fenics objects using the xdmf file

mesh, dx, ds, n, domain_meshFunction, boundary_meshFunction = file_tools.read_xdmfMesh(
file_name)

"""

# Initializes the mesh object and reads the xdmf file

mesh = Mesh()

# Initializes a mesh value collection to store mesh data

data_meshCollection = MeshValueCollection("size_t", mesh, mesh.topology(
).dim())

# Reads the mesh with domain physical groups

with XDMFFile(file_name+"_domain.xdmf") as infile:

    infile.read(mesh)

    infile.read(data_meshCollection, "domain")

ct = MeshFunction("size_t", mesh, data_meshCollection)

data_meshCollection = MeshValueCollection("size_t", mesh, mesh.topology(
).dim()-1)

# Reads the mesh with surface physical groups

with XDMFFile(file_name+"_boundary.xdmf") as infile:
   
    infile.read(data_meshCollection, "boundary")

# Converts the mesh value collections to mesh functions, for mesh value
# collections are low level and cannot be used for FEM integration and 
# other higher level operations inside FEniCs

ft = MeshFunction("size_t", mesh, data_meshCollection)

# Sets the integration differentials, with domains of integration. The 
# differentials have reserved number for each physical group:
#
# Volumetric physical groups:
# dx(1) -> fibers outside the RVE
# dx(2) -> matrix outside the RVE
# dx(3) -> RVE_fiber
# dx(4) -> RVE_matrix
#
# Surface physical groups:
# ds(5) -> bottom
# ds(6) -> front
# ds(7) -> left
# ds(8) -> back
# ds(9) -> right
# ds(10) -> top

dx = Measure("dx", domain=mesh, subdomain_data=ct)

ds = Measure("ds", domain=mesh, subdomain_data=ft)

# Sets the normal vector to the mesh's boundary

n  = FacetNormal(mesh)
"""

########################################################################
#                            Function space                            #
########################################################################

# Defines the finite element spaces for the displacement field, u, and 
# for the microrotation field, phi

U = VectorElement("CG", mesh.ufl_cell(), polynomial_degree)

V = VectorElement("CG", mesh.ufl_cell(), polynomial_degree)

# Defines the finite element spaces for the Lagrange multipliers, which 
# enforce the micro-scale kinematical constraints. The degree of the po-
# lynomial is 0 because the lagrange multiplier is constant throughout 
# the element

# Displacement constraint (impedes rigid body translation)

lagrangeMult_U = VectorElement("Real", mesh.ufl_cell(), 0) 

# Displacement gradient constraint (impedes rigid body rotation)

lagrangeMult_GradU = TensorElement("Real", mesh.ufl_cell(), 0) 

# Microrotation constraint

lagrangeMult_Phi = VectorElement("Real", mesh.ufl_cell(), 0) 

# Microrotation gradient constraint

lagrangeMult_GradPhi = TensorElement("Real", mesh.ufl_cell(), 0)

# Defines the finite element spaces for post-processing

Space_Piola_1st = TensorElement("DG", mesh.ufl_cell(), 0)

# Defines the mixed element for the monolithic solution

MicroPolarMixedElement = MixedElement([U,V, lagrangeMult_U, 
lagrangeMult_GradU, lagrangeMult_Phi, lagrangeMult_GradPhi])

monolithic_functionSpace = FunctionSpace(mesh, MicroPolarMixedElement)

# Defines the function space for the stress tensor

W = FunctionSpace(mesh, Space_Piola_1st)

########################################################################
#                          Micropolar tensors                          #
########################################################################

# Constructs the three-dimensional identity tensor 

I = Identity(3)

# Defines the curvature tensor

def curvature_tensor(phi):

    # Defines the permutation tensor for 3D euclidean space

    perm = ufl.PermutationSymbol(3)
    
    # Computes the micro-rotation tensor Rbar
     
    R_bar = t_tools.rotation_tensorEulerRodrigues(phi)
    
    # Computes the gradient of each component of R_bar to get a third 
    # order tensor, grad_Rbar

    grad_Rbar = grad(R_bar)

    # Constructs the third order curvature tensor K_third

    K_third = as_tensor(R_bar[j,i]*grad_Rbar[j,k,l], (i,k,l))

    # Contracts with the permutation tensor to form K_second

    K_second = 0.5*as_tensor(perm[i,j,k]*K_third[k,j,l], (i,l))

    return K_second

########################################################################
#                        Constitutive modelling                        #
########################################################################

# Sets the Helmholtz free energy density to a neo-hookean isotropic mo-
# del. Uses the micropolar stretch as input

def psi_NH(V_bar):

    energy = 0.5*mu*(tr(V_bar*V_bar.T)-3.0)

    return energy

# Sets the volumetric part of the Helmholtz free energy density

def psi_vol(Vbar):

    J = det(Vbar)

    energy = ((lmbda*0.25)*((J**2)-1))-((lmbda*0.5)*ln(J))-(mu*ln(J))

    return energy

# Sets the energy density related to the curvature

def psi_k(K_curvature):

    energy = ((0.5*(alpha*(tr(K_curvature)**2)+(beta*tr(K_curvature*
    K_curvature))+(gamma*tr(K_curvature*K_curvature.T)))))

    return energy

def psi_hat(Vbar):

    energy = 0.25*kappa*(tr(Vbar*Vbar.T)- tr(Vbar*Vbar)) 

    return energy

# Defines the micropolar Kirchhoff stress tensor

def Kirchhoff_Stress(u, phi):

    # Sets the identity tensor and the deformation gradient

    I = Identity(3)

    F = grad(u) + I

    # Evaluates the rotation tensor given by the micropolar rotation

    R_bar = t_tools.rotation_tensorEulerRodrigues(phi)

    # Evaluates the micropolar stretch and the jacobian
    
    V_bar = F*(R_bar.T)

    J = det(V_bar)

    # Evaluates the curvature tensor and its push-forward

    K_curvature = curvature_tensor(phi)

    k_curvatureSpatial = R_bar*K_curvature*R_bar.T

    V_barTransposed = variable(V_bar.T)

    # Evaluates the total energy density

    psi_total = (psi_NH(V_barTransposed.T)+psi_vol(V_barTransposed.T)+
    psi_hat(V_barTransposed.T)+psi_k(k_curvatureSpatial))

    # Evaluates the micropolar Kirchhoff stress

    tau = V_bar*diff(psi_total,V_barTransposed)

    # Pulls back to the reference configuration

    return J*tau*inv(def_grad(u)).T

# Definition of the Kirchhoff couple stress 

def Couple_Kirchhoff_Stress(u, phi):

    # Sets the identity tensor and the deformation gradient

    I = Identity(3)

    F = grad(u) + I

    # Evaluates the rotation tensor given by the micropolar rotation

    R_bar = t_tools.rotation_tensorEulerRodrigues(phi)

    # Evaluates the micropolar stretch and the jacobian
    
    V_bar = F*(R_bar.T)

    J = det(V_bar)

    # Evaluates the curvature tensor and its push-forward

    K_curvature = curvature_tensor(phi)

    k_curvatureSpatial = R_bar*K_curvature*R_bar.T
    
    k_curvatureSpatialTransposed = variable(k_curvatureSpatial.T)

    # Evaluates the total energy density

    psi_total = (psi_NH(V_bar)+psi_vol(V_bar)+psi_hat(V_bar)+psi_k(
    k_curvatureSpatialTransposed.T))

    # Evaluates the couple micropolar Kirchhof stress

    tau = V_bar*diff(psi_total,k_curvatureSpatialTransposed)

    # Pulls back to the reference configuration

    return J*tau*inv(def_grad(u)).T

# Definition of the deformation gradient
def def_grad(u):
    I = Identity(3)
    return I + grad(u)

########################################################################
#                       Macro quantities reading                       #
########################################################################

u_RVE_homogenized = file_tools.txt_toList(macro_displacementName)

grad_u_RVE_homogenized = file_tools.txt_toList(
macro_gradDisplacementName)

phi_RVE_homogenized = file_tools.txt_toList(macro_microrotationName)

grad_phi_RVE_homogenized = file_tools.txt_toList(
macro_gradMicrorotationName)

########################################################################
#                         Material properties                          #
########################################################################

# Converts the dictionary of shear modulus to a discontinuous Galerkin
# function space

mu = fun_tools.physical_groupToDGSpace(mu_materials, mesh, 
domain_meshFunction)

# Defines other material properties derived from the shear modulus (the
# mu Lamé parameter)

kappa = -0.15*mu

lmbda = 63.84-(2*mu/3)

# Calculates the volume of the RVE

v_tot = assemble(1*dx)

########################################################################
#                           Variational forms                          #
########################################################################

# Defines the trial and test functions

dsol = TrialFunction(monolithic_functionSpace) 

v = TestFunction(monolithic_functionSpace)

# Splits the test functions into their respective fields

(v_displacement, v_microrotation, v_lagrangeMultDisplacement, 
v_lagrangeMultGradDisp, v_lagrangeMultMicrorotation, 
v_lagrangeMultGradRotation) = split(v)

# Creates the function for the updated solution, i.e. the vector of pa-
# rameters

sol_new = Function(monolithic_functionSpace)

# Splits the updated solution into their respective fields

(new_displacement, new_microrotation, new_lagrangeDisplacement, 
new_lagrangeGradDisp, new_lagrangeRotation, new_lagrangeGradRotation
) = split(sol_new)

# Turns the stresses into objects for easy handling

stress = Kirchhoff_Stress(new_displacement,new_microrotation)

couple_stress = Couple_Kirchhoff_Stress(new_displacement,
new_microrotation)

# Defines the variational forms 

def BVP(time_step):

    # Gets the current time step from the file of the macro displacement

    t_current = u_RVE_homogenized[time_step][0]

    # Converts the macro quantities to vectors and tensors

    macro_displacement = as_vector(u_RVE_homogenized[time_step][1])

    macro_rotation = as_vector(phi_RVE_homogenized[time_step][1])

    macro_gradDisp = as_tensor(grad_u_RVE_homogenized[time_step][1])

    macro_gradRotation = as_tensor(grad_phi_RVE_homogenized[time_step][1])

    # Defines the variational forms

    Int_du = inner(stress, grad(v_displacement))*dx 

    Int_dphi = ((inner(couple_stress,grad(v_microrotation))*dx)-(inner(
    stress*def_grad(new_displacement).T,t_tools.skew_2OrderTensor(
    v_microrotation))*dx)) 

    Int_dlambu = ((1/v_tot)*((dot(v_lagrangeMultDisplacement,
    macro_displacement)*dx)-(dot(v_lagrangeMultDisplacement,
    new_displacement)*dx)))

    Int_dlambphi = ((1/v_tot)*((dot(v_lagrangeMultMicrorotation,
    macro_rotation)*dx)-(dot(v_lagrangeMultMicrorotation,
    new_microrotation)*dx)))

    Int_dlambgrau = ((1/v_tot)*((inner(v_lagrangeMultGradDisp, 
    macro_gradDisp)*dx)-(inner(v_lagrangeMultGradDisp, grad(
    new_displacement))*dx)))

    Int_dlambgradphi = ((1/v_tot)*((inner(v_lagrangeMultGradRotation, 
    macro_gradRotation)*dx)-(inner(v_lagrangeMultGradRotation,grad(
    new_microrotation))*dx))) 

    l_du = (-(1/v_tot)*((dot(new_lagrangeDisplacement, v_displacement)*
    dx)+(inner(new_lagrangeGradDisp, grad(v_displacement))*dx)))

    l_dphi = (-(1/v_tot)*((dot(new_lagrangeRotation,v_microrotation)*dx)
    +(inner(new_lagrangeGradRotation, grad(v_microrotation))*dx))) 

    # Sums all variational forms to evaluate the residual

    residual = (((1/v_tot)*(Int_du+Int_dphi))+Int_dlambu+Int_dlambgrau+
    Int_dlambphi+Int_dlambgradphi+l_du+l_dphi)

    # Evaluates the jacobian of the residual

    derivative_residual = derivative(residual, sol_new, dsol)

    # Defines the nonlinear variational problem and solver parameters

    Res = NonlinearVariationalProblem(residual, sol_new, None, J=
    derivative_residual)

    solver = NonlinearVariationalSolver(Res)

    solver.parameters["nonlinear_solver"] = "newton"

    solver.parameters["newton_solver"]["linear_solver"] = linear_solver

    solver.parameters["newton_solver"]["relative_tolerance"] = (
    relative_tolerance)

    solver.parameters["newton_solver"]["maximum_iterations"] = (
    maximum_iterations)
    
    solver.parameters["newton_solver"]["preconditioner"] = (
    preconditioner)

    solver.parameters['newton_solver']['krylov_solver']['absolute_tole'+
    'rance'] = krylov_absoluteTolerance

    solver.parameters['newton_solver']['krylov_solver']['relative_tole'+
    'rance'] = krylov_relativeTolerance

    solver.parameters['newton_solver']['krylov_solver']['maximum_itera'+
    'tions'] = krylov_maximumIterations

    solver.parameters['newton_solver']['krylov_solver']['monitor_conve'+
    'rgence'] = krylov_monitorConvergence

    return solver, t_current

########################################################################
#                         Files initialization                         #
########################################################################

# Initializes the XMDF files

microrotation_file = XDMFFile(MPI.COMM_WORLD, "./ResultsDir/phi_micro."+
"xdmf")

displacement_file = XDMFFile(MPI.COMM_WORLD, "./ResultsDir/u_micro.xdm"+
"f")

first_piolaFile = XDMFFile(MPI.COMM_WORLD, "./ResultsDir/PK1_micro.xdm"+
"f")

couple_firstPiolaFile = XDMFFile(MPI.COMM_WORLD, "./ResultsDir/PK1_cou"+
"ple_micro.xdmf")

# Opens and makes the headers for the homogenized quantities and for the
# Lagrange multipliers

file_LagrangeMultiplier_Stress = open("./ResultsDir/LagrangeMultiplier"+
"_Stress.txt", "w")

file_LagrangeMultiplier_Stress.write("Time\tP11\tP12\tP13\tP21\tP22\tP"+
"23\tP31\tP32\tP33\n")  

file_LagrangeMultiplier_Couple_Stress = open("./ResultsDir/LagrangeMul"+
"tiplier_Couple_Stress.txt", "w")

file_LagrangeMultiplier_Couple_Stress.write("Time\tPc11\tPc12\tPc13\tP"+
"c21\tPc22\tPc23\tPc31\tPc32\tPc33\n")  

file_Homogenized_Stress = open("./ResultsDir/Homogenized_Stress.txt", 
"w")

file_Homogenized_Stress.write("Time\tP11\tP12\tP13\tP21\tP22\tP23\tP31"+
"\tP32\tP33\n")  

file_Homogenized_Couple_Stress = open("./ResultsDir/Homogenized_Couple"+
"d_Stress.txt", "w")

file_Homogenized_Couple_Stress.write("Time\tPc11\tPc12\tPc13\tPc21\tPc"+
"22\tPc23\tPc31\tPc32\tPc33\n")  

########################################################################
#                   Solution and pseudotime stepping                   #
########################################################################

# Initializes the pseudo time counter

time_counter = 0

# Initializes functions for the first Piola-Kirchhoff stress tensors, so
# that they have a finite element interpolation

PK1 = Function(W, name='PK1')

PK1_couple = Function(W, name='PK1_couple')

# Iterates through the pseudotime stepping

while True:

    # Retrieves the boundary value problem information

    solver, t_current = BVP(time_counter)

    print("###########################################################"+
    "#############\n#                 Incremental step: "+str(
    time_counter+1)+"; current time: "+str(t_current)+"               "+
    "#\n##############################################################"+
    "##########\n")

    # Solves the nonlinear variational problem

    solver.solve()

    # Splits the solution in their corresponding fields

    u, phi, ulag, gradulag, philag, gradphilag = sol_new.split(deepcopy=
    True)

    u.rename("Micro-scale Displacement", "Multiscale")
    phi.rename("Micro-scale Micro-rotation", "Multiscale")
    
    # Saves the values of the Lagrange multiplier relative to the gra-
    # dient of the displacement

    gradulag_value_at_t = gradulag.vector().get_local()

    # Extracts its components
    T11, T12, T13, T21, T22, T23, T31, T32, T33 = gradulag_value_at_t[0:9]

    # Writes time and tensor components to file

    file_LagrangeMultiplier_Stress.write(f"{t_current:.3f}\t{T11:.6f}\t{T12:.6f}\t{T13:.6f}\t{T21:.6f}\t{T22:.6f}\t{T23:.6f}\t{T31:.6f}\t{T32:.6f}\t{T33:.6f}\n")

    # Extract the values of gradphilag and save them to the file
    gradphilag_value_at_t = gradphilag.vector().get_local()
    # Extract components
    Tc11, Tc12, Tc13, Tc21, Tc22, Tc23, Tc31, Tc32, Tc33 = gradphilag_value_at_t[0:9]  # Adjust based on tensor shape
    # Write time and tensor components to file
    file_LagrangeMultiplier_Couple_Stress.write(f"{t_current:.3f}\t{Tc11:.6f}\t{Tc12:.6f}\t{Tc13:.6f}\t{Tc21:.6f}\t{Tc22:.6f}\t{Tc23:.6f}\t{Tc31:.6f}\t{Tc32:.6f}\t{Tc33:.6f}\n")

    displacement_file.write(u, t_current)
    microrotation_file.write(phi, t_current)

    PK1 = project(Kirchhoff_Stress(u,phi), W)
    PK1.rename("Micro-scale PK1 Stress", "Multiscale")
    first_piolaFile.write(PK1, t_current)
    PK1_couple = project(Couple_Kirchhoff_Stress(u,phi), W)
    PK1_couple.rename("Micro-scale PK1 Couple Stress", "Multiscale")
    couple_firstPiolaFile.write(PK1_couple, t_current)

    Homogenized_Piola11 = assemble(Kirchhoff_Stress(u,phi)[0,0]*dx)/v_tot
    Homogenized_Piola12 = assemble(Kirchhoff_Stress(u,phi)[0,1]*dx)/v_tot
    Homogenized_Piola13 = assemble(Kirchhoff_Stress(u,phi)[0,2]*dx)/v_tot
    Homogenized_Piola21 = assemble(Kirchhoff_Stress(u,phi)[1,0]*dx)/v_tot
    Homogenized_Piola22 = assemble(Kirchhoff_Stress(u,phi)[1,1]*dx)/v_tot
    Homogenized_Piola23 = assemble(Kirchhoff_Stress(u,phi)[1,2]*dx)/v_tot
    Homogenized_Piola31 = assemble(Kirchhoff_Stress(u,phi)[2,0]*dx)/v_tot
    Homogenized_Piola32 = assemble(Kirchhoff_Stress(u,phi)[2,1]*dx)/v_tot
    Homogenized_Piola33 = assemble(Kirchhoff_Stress(u,phi)[2,2]*dx)/v_tot

    # Extract the values of gradulag and save them to the file
    # Extract components
    T11 = Homogenized_Piola11 # Adjust based on tensor shape
    T12 = Homogenized_Piola12 # Adjust based on tensor shape
    T13 = Homogenized_Piola13 # Adjust based on tensor shape
    T21 = Homogenized_Piola21 # Adjust based on tensor shape
    T22 = Homogenized_Piola22 # Adjust based on tensor shape
    T23 = Homogenized_Piola23 # Adjust based on tensor shape
    T31 = Homogenized_Piola31 # Adjust based on tensor shape
    T32 = Homogenized_Piola32 # Adjust based on tensor shape
    T33 = Homogenized_Piola33 # Adjust based on tensor shape

    # Write time and tensor components to file
    file_Homogenized_Stress.write(f"{t_current:.3f}\t{T11:.6f}\t{T12:.6f}\t{T13:.6f}\t{T21:.6f}\t{T22:.6f}\t{T23:.6f}\t{T31:.6f}\t{T32:.6f}\t{T33:.6f}\n")

    time_counter += 1

    if time_counter>=maximum_loadingSteps:

        print("The maximum number of loading steps,",
        maximum_loadingSteps, "has just been reached. Stops the simula"+
        "tion immediatly\n")

        break

file_LagrangeMultiplier_Stress.close()
file_LagrangeMultiplier_Couple_Stress.close()

print ('Simulation Completed')
# Turn off interactive mode when done