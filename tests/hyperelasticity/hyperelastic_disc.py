# Hyperelasticity

from fenics import *

# Form compiler options
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True

##################################################
#                Pseudotime control              #
##################################################

# Defines materials properties

E_young = 200E6

nu = 0.3

mu = E_young/(2*(1+nu))

lmbda = (E_young*nu)/((1+nu)*(1-2*nu))

t = 0.0

t_final = 1.0

n_steps = 10

delta_t = t_final/n_steps

# Defines a mesh

z_max = 0.1

mesh = BoxMesh(Point(0.0,0.0,0.0),Point(0.01,0.01,z_max),2,2,15)

# Defines a vectorial function space

V = VectorFunctionSpace(mesh, "P", 2)

# Defines a trial function

Delta_u = TrialFunction(V)

u = Function(V)

# Defines a test function

du = TestFunction(V)

# Defines Dirichlet boundary conditions

#uz_max = 0.0*z_max

below = CompiledSubDomain('near(x[2],z_boundary) && on_boundary', z_boundary=0.0)

above = CompiledSubDomain('near(x[2],z_boundary) && on_boundary', z_boundary=z_max)

#expression_above = Expression(("(t/t_final)*uz_max"),
#t_final=t_final, uz_max=uz_max, t=0.0, degree=1)

expression_below = Constant(("0.0", "0.0", "0.0"))

bc_dirichlet = [DirichletBC(V, expression_below, below)]#,DirichletBC(V.sub(2), expression_above, above)]

# Defines Neumann boundary conditions

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1,0)

above.mark(boundaries,1)

ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Defines the traction vector

T_max = 0.0#1E6

T = Expression(("0.0", "(t/t_final)*T_max", "0.0"), t_final=t_final,
t=0.0, T_max=T_max, degree=0) 

# Defines the body forces vector

B = Constant(("0.0","0.0","0.0"))

###################################################
#                    Tensorial                    #
###################################################

I = Identity(3)

u.interpolate(Constant((0,0,0)))

# Defines the deformation gradient

F = grad(u)+I 

# Defines the Euler-Lagrange strain tensor

C = variable((F.T)*F) 

E = 0.5*(C-I)

# Defines the Saint Venant energy function

#psi = (mu*inner(E,E))+(0.5*lmbda*(tr(E)**2))

I1 = tr(C)

J = sqrt(det(C))

psi = mu/2*(I1 - 3) + lmbda*((0.5*pow((J-1), 2))-ln(J))# = 0.5*mu*(I1-3)+0.5*lmbda*(pow((J-1),2))

# Defines the second Piola-Kirchhof stress tensor

S = 2*diff(psi, C)#(lmbda*tr(E)*I)+(2*mu*E)#2*diff(psi, C)

# Defines the first Piola-Kirchhof stress tensor

P = F*S

dF = grad(du)

###################################################
#                  Variational                    #
# #################################################

# Defines the bilinear form

a = inner(P,dF)*dx

# Defines the linear form

L = (dot(B,du)*dx)+(dot(T,du)*ds(1))

residue = a-L

Jacobian = derivative(residue, u, Delta_u)

problem = NonlinearVariationalProblem(residue, u, bc_dirichlet, J=Jacobian)

solution = NonlinearVariationalSolver(problem)

solution.parameters["nonlinear_solver"] = "newton"

solution.parameters["newton_solver"]["absolute_tolerance"] = 1E-8

solution.parameters["newton_solver"]["relative_tolerance"] = 1E-8

solution.parameters["newton_solver"]["maximum_iterations"] = 50

solution.parameters["newton_solver"]["linear_solver"] = "mumps"

# Sweeps through the pseudotime steps

reaction_z = 0.0

# Creates the displacement file

file = File("displacement.pvd")

for i in range(n_steps):

    print("\nRuns pseudotime", i, ":", t, "\n")

    T.t = t

    #expression_above.t = t

    # Solves the system

    solution.solve()

    N = as_vector([0.0,0.0,1.0])

    reaction_z = assemble((dot(N,P*N))*ds)

    print("\nThe reaction in the z direction is:",
    reaction_z, "\n")

    file << u

    # Updates time and the traction vector

    t += delta_t