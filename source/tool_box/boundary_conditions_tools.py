# Routine to store methods to neatly apply boundary conditions

from dolfin import *

########################################################################
#               Homogeneous Dirichlet boundary conditions              #
########################################################################

# Defines a function to apply homogeneous boundary conditions to a vec-
# tor field over all the 3 directions. The field can be mixed, i.e. the-
# re can be multiple field within the same solution. The number of 
# fields must be provided because, if one programatically tries to count
# the number of subfields in a function space and the function space is
# a vector function space for instance, it will count the dimensionality
# of the vector as subfields, which are not physical fields alone per
# se. If sub_fields is an empty list, all DOFs are constrained

def fixed_supportDirichletBC(field_functionSpace, boundary_meshFunction, 
boundary_physicalGroups=0, n_fields=1, sub_fieldsToApplyBC=[]):

    # Initializes a list of boundary conditions objects

    boundary_conditions = []

    # If the physical groups variable is null, returns the empty list of
    # boundary conditions

    if boundary_physicalGroups==0:

        return boundary_conditions

    # Verifies whether the boundary physical groups is a list or not

    if isinstance(boundary_physicalGroups, list):

        # Iterates through the regions

        for physical_group in boundary_physicalGroups:

            # Verifies if there is only one field

            if n_fields==1:

                # Adds this particular boundary condition to the lot

                boundary_conditions.append(DirichletBC(
                field_functionSpace, Constant((0.0, 0.0, 0.0)), 
                boundary_meshFunction, physical_group))

            # Otherwise, iterates through the fields

            else:

                for field in range(n_fields):

                    # Verifies if this field is to be constrained or not

                    if (field in sub_fieldsToApplyBC or (
                    sub_fieldsToApplyBC==[])):

                        # Adds this particular boundary condition to the 
                        # lot

                        boundary_conditions.append(DirichletBC(
                        field_functionSpace.sub(field), Constant((0.0, 0.0, 
                        0.0)), boundary_meshFunction, physical_group))

    # Otherwise, if there is only one physical group to apply boundary 
    # conditions

    else:

        # Verifies if there is only one field

        if n_fields==1:

            # Adds this particular boundary condition to the lot

            boundary_conditions.append(DirichletBC(field_functionSpace,
            Constant((0.0, 0.0, 0.0)), boundary_meshFunction, 
            boundary_physicalGroups))

        # Otherwise, iterates through the fields

        else:

            for field in range(n_fields):

                # Verifies if this field is to be constrained or not

                if (field in sub_fieldsToApplyBC or (
                sub_fieldsToApplyBC==[])):

                    # Adds this particular boundary condition to the lot

                    boundary_conditions.append(DirichletBC(
                    field_functionSpace.sub(field), Constant((0.0, 0.0, 
                    0.0)), boundary_meshFunction, 
                    boundary_physicalGroups))

    # Returns the boundary conditions list

    return boundary_conditions

# Defines a function to apply a simple support. If sub_fields is an empty
# list, all DOFs are constrained

def simple_supportDirichletBC(field_functionSpace, boundary_meshFunction, 
boundary_physicalGroups, list_constrainedDOFs, n_fields=1, 
sub_fieldsToApplyBC=[]):

    # Initializes a list of boundary conditions objects

    boundary_conditions = []

    # Verifies whether the boundary physical groups is a list or not

    if isinstance(boundary_physicalGroups, list):

        # Iterates through the regions

        for physical_group in boundary_physicalGroups:

            # Verifies if there is only one field

            if n_fields==1:

                # Iterates through the DOFs to be constrained

                for i in range(field_functionSpace.num_sub_spaces()):
                    
                    if i in list_constrainedDOFs:

                        # Adds this particular boundary condition to the 
                        # lot

                        boundary_conditions.append(DirichletBC(
                        field_functionSpace.sub(i), Constant(0.0), 
                        boundary_meshFunction, physical_group))

            # Otherwise, iterates through the fields

            else:

                for field in range(n_fields):

                    # Verifies if this field is to be constrained or not

                    if (field in sub_fieldsToApplyBC or (
                    sub_fieldsToApplyBC==[])):
                        
                        # Iterates through the DOFs to be constrained

                        for i in range(field_functionSpace.sub(field
                        ).num_sub_spaces()):
                            
                            if i in list_constrainedDOFs:

                                # Adds this particular boundary condition 
                                # to the lot

                                boundary_conditions.append(DirichletBC(
                                field_functionSpace.sub(field).sub(i), 
                                Constant(0.0), boundary_meshFunction, 
                                physical_group))

    # Otherwise, if there is only one physical group to apply boundary 
    # conditions

    else:

        # Verifies if there is only one field

        if n_fields==1:

            # Iterates through the DOFs to be constrained

            for i in range(field_functionSpace.num_sub_spaces()):
                
                if i in list_constrainedDOFs:

                    # Adds this particular boundary condition to the lot

                    boundary_conditions.append(DirichletBC(
                    field_functionSpace.sub(i), Constant(0.0), 
                    boundary_meshFunction, boundary_physicalGroups))

        # Otherwise, iterates through the fields

        else:

            for field in range(n_fields):

                # Verifies if this field is to be constrained or not

                if (field in sub_fieldsToApplyBC or (
                sub_fieldsToApplyBC==[])):
                    
                    # Iterates through the DOFs to be constrained

                    for i in range(field_functionSpace.sub(field
                    ).num_sub_spaces()):
                        
                        if i in list_constrainedDOFs:

                            # Adds this particular boundary condition to
                            # the lot

                            boundary_conditions.append(DirichletBC(
                            field_functionSpace.sub(field).sub(i), 
                            Constant(0.0), boundary_meshFunction, 
                            boundary_physicalGroups))

    # Returns the boundary conditions list

    return boundary_conditions