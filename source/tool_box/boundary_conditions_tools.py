# Routine to store methods to neatly apply boundary conditions

from dolfin import *

########################################################################
#               Homogeneous Dirichlet boundary conditions              #
########################################################################

# Defines a function to apply homogeneous boundary conditions to a vec-
# tor field over all the 3 directions. The field can be mixed, i.e. the-
# re can be multiple field within the same solution

def cantilever_DirichletBC(field_function, supported_dofsList, 
boundary_meshFunction, boundary_physicalGroups, n_fields=1):

    # Initializes a list of boundary conditions objects

    boundary_conditions = []

    # Verifies whether the boundary physical groups is a list or not

    if isinstance(boundary_physicalGroups, list):

        #

    bc1 = DirichletBC(monolithic_functionSpace.sub(0), Constant((0.0, 0.0, 0.0)), boundary_meshFunction, 5)
    bc2 = DirichletBC(monolithic_functionSpace.sub(1), Constant((0.0, 0.0, 0.0)), boundary_meshFunction, 5)