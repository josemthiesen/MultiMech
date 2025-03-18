# Routine to compute the micropolar macroscale problem

from dolfin import *

from mpi4py import MPI

import ufl_legacy as ufl

import numpy as np

import matplotlib.pyplot as plt

from mshr import *

#import periodic_structure as mesher

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.homogenization_tools as homogenization_tools

import source.tool_box.file_handling_tools as file_tools

import source.tool_box.functional_tools as functional_tools

import source.tool_box.boundary_conditions_tools as BCs_tools

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

gamma = 1e-12#5.22e-8 #5.22e-8 #((1e-5)**2)*2*mu # 43.265e-6

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

file_name = "test_meshes//macro_mesh"

# Defines a flag to generate a new mesh or not

flag_newMesh = False

# Defines the physical groups that belongs to the RVE

volume_physGroupsRVE = [3, 4]

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

# Reads the mesh and constructs some fenics objects using the xdmf file

(mesh, dx, ds, n, domain_meshCollection, domain_meshFunction, 
boundary_meshCollection, boundary_meshFunction) = mesh_tools.read_xdmfMesh(
file_name)

# Defines the finite element spaces for the displacement field, u, and 
# for the microrotation field, phi

U = VectorElement("CG", mesh.ufl_cell(), polynomial_degree)

V = VectorElement("CG", mesh.ufl_cell(), polynomial_degree)

# Define the mixed element for the monolithic solution

micropolar_mixedElement = MixedElement([U,V])
 
monolithic_functionSpace = FunctionSpace(mesh, micropolar_mixedElement)

########################################################################
#                         RVE post-processing                          #
########################################################################

# Creates the RVE submesh and automatically constructs the function spa-
# ces and the DOFs mapping between meshes

(RVE_submesh, domain_meshFunction, UV_submesh, RVE_meshMapping, 
parent_meshMapping, sol_RVE, RVE_toParentCellMap) = mesh_tools.create_mesh(
mesh, domain_meshFunction, volume_physGroupsRVE, 
monolithic_functionSpace, mixed_element=micropolar_mixedElement)

########################################################################
#                         Material properties                          #
########################################################################

# Converts the dictionary of shear modulus to a discontinuous Galerkin
# function space

mu = functional_tools.physical_groupToDGSpace(mu_materials, mesh, 
domain_meshFunction)

# Defines other material properties derived from the shear modulus (the
# mu Lamé parameter)

kappa = 0#-0.15*mu

lmbda = 63.84-(2*mu/3)

########################################################################
#                          Boundary conditions                         #
########################################################################

# Defines a load expression

load = Expression("t", t = 0.0, degree = 1)

# Defines the boundary conditions

bc = BCs_tools.fixed_supportDirichletBC(monolithic_functionSpace,)

bc1 = DirichletBC(monolithic_functionSpace.sub(0), Constant((0.0, 0.0, 0.0)), boundary_meshFunction, 5)
bc2 = DirichletBC(monolithic_functionSpace.sub(1), Constant((0.0, 0.0, 0.0)), boundary_meshFunction, 5)

bc = [bc1, bc2]

traction_boundary = as_vector([load,0.0,0.0])

I = Identity(3)

def safe_sqrt(a):
    return sqrt(a + 1.0e-12)


def phi_angle(phi):
    return safe_sqrt(dot(phi, phi))

def W_matrix(phi):
   
   perm = ufl.PermutationSymbol(3)
   
   W = as_tensor(perm[j,i,k]*phi[k], (i,j))

   return W

def MicroRotation_Tensor(phi):

    rotation_angle = phi_angle(phi)

    W_phi = W_matrix(phi)

    R_bar = ufl.cos(rotation_angle)*I + (ufl.sin(rotation_angle)/rotation_angle)*W_phi \
                + ((1-ufl.cos(rotation_angle))/(rotation_angle**2))*ufl.outer(phi,phi)

    return R_bar

def curvature_tensor(phi):

    # Define the permutation tensor for 3D
    perm = ufl.PermutationSymbol(3)
    
    # Compute the micro-rotation tensor `Rbar` 
    Rbar = MicroRotation_Tensor(phi)
    
    # Compute the gradient of each component of Rbar to get a third-order tensor `gradRbar`
    # Here we assume that Rbar has dimensions (3, 3), and we calculate the gradient manually
    gradRbar = grad(Rbar)

    # Construct the third-order curvature tensor Kthird
    Kthird = as_tensor(Rbar[j, i] * gradRbar[j,k,l], (i, k, l))

    # Contract with the permutation tensor to form Ksecond
    Ksecond = (1/2)*as_tensor(perm[i, j, k] * Kthird[k, j, l], (i, l))

    return Ksecond

def psi_vol(Vbar):

    J = det(Vbar)

    energy = ((lmbda/4)*((J**2) - 1)) - ((lmbda/2)*(ln(J))) - (mu*ln(J))

    return energy

def psi_NH(Vbar):

    energy = 0.5*mu*(tr(Vbar*Vbar.T)-3)
    return energy

def psi_hat(Vbar):

    energy = 0.25*kappa*(tr(Vbar*Vbar.T)- tr(Vbar*Vbar)) 

    return energy

def psi_k(K_curvature):

    energy = 0.5*(alpha*(tr(K_curvature)**2) + beta*tr(K_curvature*K_curvature) + gamma*tr(K_curvature*K_curvature.T))

    return energy


def Kirchhoff_Stress(u, phi):

    I = Identity(3)
    F = grad(u) + I

    R_bar = MicroRotation_Tensor(phi)
    
    V_bar = F*(R_bar.T)
    J = det(V_bar)

    #K_curvature = curvature_tensor(phi)
    #k_curvature_spatial = R_bar*K_curvature*R_bar.T

    #V_bar_trans = variable(V_bar.T)
    #k_curvature_spatial = variable(k_curvature_spatial.T)

    #psi_total = psi_NH(V_bar_trans.T) + psi_vol(V_bar_trans.T) +  psi_hat(V_bar_trans.T) + psi_k(k_curvature_spatial)

    #tau = J*V_bar*diff(psi_total, V_bar_trans)

    tau = ((lmbda/2)*((J*J)-1)*I) + (mu*((V_bar*V_bar.T)-I)) + ((kappa/2)*((V_bar*V_bar.T)-(V_bar*V_bar)))

    return tau

def Couple_Kirchhoff_Stress(u, phi):

    I = Identity(3)
    F = grad(u) + I
    
    R_bar = MicroRotation_Tensor(phi)
    
    V_bar = F*(R_bar.T)
    J = det(V_bar)

    K_curvature = curvature_tensor(phi)
    k_curvature_spatial = R_bar*K_curvature*R_bar.T

    #V_bar = variable(V_bar.T)
    #k_curvature_spatial_trans = variable(k_curvature_spatial.T)

    #psi_total = psi_NH(V_bar) + psi_vol(V_bar) +  psi_hat(V_bar) + psi_k(k_curvature_spatial_trans.T)

    #tau = J*V_bar*diff(psi_total, k_curvature_spatial_trans)

    tau = V_bar*((alpha*tr(k_curvature_spatial)*I)+(beta*k_curvature_spatial) + (gamma*k_curvature_spatial.T))

    return tau

def def_grad(u):
    I = Identity(3)
    return I + grad(u)

# The function that multiplied the PDE/residual is called test function. The unknown function u to be approximated is referred to as trial function.

dsol, v = TrialFunction(monolithic_functionSpace), TestFunction(monolithic_functionSpace)
vu, vphi = split(v)

sol_new = Function(monolithic_functionSpace)

u_new, phi_new = split(sol_new)

invgrad = inv(def_grad(u_new)).T

stress = Kirchhoff_Stress(u_new,phi_new)
couple_stress = Couple_Kirchhoff_Stress(u_new,phi_new)

Int_du = inner(stress*invgrad, grad(vu))*dx 

Int_dphi = inner(couple_stress*invgrad, grad(vphi))*dx - inner(stress, W_matrix(vphi))*dx 

l_du = dot(traction_boundary,vu)*ds(9)

#l_dphi = dot(moment,vphi)*ds(10)

F_Int = Int_du + Int_dphi - l_du #- l_dphi

Gain = derivative(F_Int , sol_new, dsol)
Res = NonlinearVariationalProblem(F_Int, sol_new, bc, J=Gain)
solver = NonlinearVariationalSolver(Res)
solver.parameters['nonlinear_solver'] = 'newton'
solver.parameters['newton_solver']['relative_tolerance'] = 1e-10
solver.parameters['newton_solver']['maximum_iterations'] = 10
solver.parameters['newton_solver']['linear_solver'] = 'mumps'
solver.parameters['newton_solver']['absolute_tolerance'] = 1e-10

conc_u = File("./ResultsDir/u.pvd")
conc_phi = File("./ResultsDir/phi.pvd")

conc_u_RVE = File("./ResultsDir/u_RVE.pvd")
conc_phi_RVE = File("./ResultsDir/phi_RVE.pvd")

iter = 0
t = 0.0
Tf = 1.0
u_r = 0.05
deltaT = 0.1

u_RVE_homogenized = []
grad_u_RVE_homogenized = []

phi_RVE_homogenized = []
grad_phi_RVE_homogenized =[]

time_counter = 0

while t < Tf:

    print("\n-------------------------------------------------------------------------")
    print("-              Incremental step:", time_counter+1, "; current time:", t, "                 -")
    print("-------------------------------------------------------------------------\n")

    solver.solve()

    u,phi = sol_new.split(deepcopy=True)

    u.rename("DNS Displacement", "DNS")
    phi.rename("DNS Micro-rotation", "DNS")

    u_RVE_homogenized.append([t,homogenization_tools.get_homogenized(u, dx).tolist()])
    grad_u_RVE_homogenized.append([t,homogenization_tools.get_homogenized_gradient(grad(u), dx).tolist()])
    phi_RVE_homogenized.append([t,homogenization_tools.get_homogenized(phi, dx).tolist()])
    grad_phi_RVE_homogenized.append([t,homogenization_tools.get_homogenized_gradient(grad(phi), dx).tolist()])

    """

    # Iterates through the elements of the submesh

    for element in cells(RVE_submesh):

        # Gets the index of the element in the RVE submesh

        submesh_index = element.index()

        # Gets the index of the element in the parent mesh

        parent_index = RVE_toParentCellMap[submesh_index]

        # Translates the values of the solution using the DOFs mapping

        sol_RVE.vector()[U_rveDofMap.cell_dofs(submesh_index)] = sol_new.vector()[U_parentDofMap.cell_dofs(parent_index)] 
        sol_RVE.vector()[V_rveDofMap.cell_dofs(submesh_index)] = sol_new.vector()[V_parentDofMap.cell_dofs(parent_index)] 
    """

    sol_RVE = mesh_tools.field_parentToSubmesh(RVE_submesh, sol_RVE, 
    sol_new, RVE_toParentCellMap, RVE_meshMapping, parent_meshMapping)

    u_RVE, phi_RVE = sol_RVE.split(deepcopy=True)

    u_RVE.rename("DNS-Subdomain Displacement", "DNS-Subdomain ")
    phi_RVE.rename("DNS-Subdomain Micro-rotation", "DNS-Subdomain ")

    conc_u << u
    conc_phi << phi
    conc_u_RVE << u_RVE
    conc_phi_RVE << phi_RVE

    t += deltaT
    load.t = t*u_r
    time_counter += 1

    file_tools.list_toTxt(u_RVE_homogenized, "u_RVE_homogenized")
    file_tools.list_toTxt(grad_u_RVE_homogenized, "grad_u_RVE_homogenized")
    file_tools.list_toTxt(phi_RVE_homogenized, "phi_RVE_homogenized")
    file_tools.list_toTxt(grad_phi_RVE_homogenized, "grad_phi_RVE_homogenized")

print ('Simulation Completed')
