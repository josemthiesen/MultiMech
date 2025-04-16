# Routine to store methods to neatly apply boundary conditions

from dolfin import *

########################################################################
#              Heterogeneous Dirichlet boundary conditions             #
########################################################################

# Defines a function to apply prescribed boundary conditions to a pro-
# blem using an Expression as load prescriber

def prescribed_DirichletBC(prescribed_conditionsDict, 
field_functionSpace, boundary_meshFunction, boundary_conditions=[], 
boundary_physGroupsNamesToTags=dict(), verbose=False):
    
    # Tests if the element is mixed and gets the number of fields

    n_fields = 1

    element_type = field_functionSpace.ufl_element().family()

    if element_type=='Mixed':

        # Gets the number of fields

        n_fields = field_functionSpace.ufl_element().num_sub_elements()

    # Iterates through the keys and values of the dictionary of prescri-
    # bed Dirichlet boundary conditions. The key is a physical group or
    # a tuple of physical groups, whereas the values are lists of: num-
    # ber of the field if there is more than one; degrees of freedom to
    # be prescribed; and the load expression

    for physical_groups, load_info in prescribed_conditionsDict.items():

        # Verifies whether the boundary physical groups is a list or not

        if isinstance(physical_groups, tuple):

            # Iterates through the regions

            for physical_group in physical_groups:

                # Verifies if the physical group is a string

                physical_group = verify_stringPhysicalGroup(physical_group, 
                boundary_physGroupsNamesToTags)

                # Verifies if there is only one field

                if n_fields==1:

                    # Verifies if no DOFs are especified

                    if len(load_info)==1:

                        # Adds this particular boundary condition to the 
                        # lot

                        boundary_conditions.append(DirichletBC(
                        field_functionSpace, load_info[0], 
                        boundary_meshFunction, physical_group))

                    # Otherwise, checks if there is a list of DOFs to be
                    # prescribed

                    else:

                        if isinstance(load_info[0], list):

                            # Iterates through the prescribed DOFs

                            for dof_number in load_info[0]:

                                # Adds this particular boundary condition 
                                # to the lot

                                boundary_conditions.append(DirichletBC(
                                field_functionSpace.sub(dof_number), 
                                load_info[1], boundary_meshFunction, 
                                physical_group))

                        else:

                            # Adds this particular boundary condition to 
                            # the lot

                            boundary_conditions.append(DirichletBC(
                            field_functionSpace.sub(load_info[0]), 
                            load_info[1], boundary_meshFunction, 
                            physical_group))

                # Otherwise, iterates through the fields

                else:

                    for field in range(n_fields):

                        # Verifies if this field is to be constrained or 
                        # not

                        if field==load_info[0]:

                            # Checks whether a set of degrees of freedom
                            # of the field is to be prescribed

                            if len(load_info)>2:

                                # Verifies if there is a list of DOFs

                                if isinstance(load_info[1], list):

                                    # Iterates through the prescribed
                                    # DOFs

                                    for dof_number in load_info[1]:

                                        # Adds this particular boundary
                                        # condition to the lot

                                        boundary_conditions.append(DirichletBC(
                                        field_functionSpace.sub(field).sub(
                                        dof_number), load_info[2], 
                                        boundary_meshFunction, 
                                        physical_group))

                                else:

                                    # Adds this particular boundary con-
                                    # dition to the lot

                                    boundary_conditions.append(DirichletBC(
                                    field_functionSpace.sub(field).sub(
                                    load_info[1]), load_info[2], 
                                    boundary_meshFunction, 
                                    physical_group))

                            else:

                                # Adds this particular boundary condi-
                                # tion to the lot

                                boundary_conditions.append(DirichletBC(
                                field_functionSpace.sub(field), 
                                load_info[1], boundary_meshFunction, 
                                physical_group))

        # Otherwise, if there is only one physical group to apply boundary 
        # conditions

        else:

            # Verifies if the physical group is a string

            physical_group = verify_stringPhysicalGroup(physical_groups, 
            boundary_physGroupsNamesToTags)

            # Verifies if there is only one field

            if n_fields==1:

                # Verifies if no DOFs are especified

                if len(load_info)==1:

                    # Adds this particular boundary condition to the lot

                    boundary_conditions.append(DirichletBC(
                    field_functionSpace, load_info[0], 
                    boundary_meshFunction, physical_group))

                # Otherwise, checks if there is a list of DOFs to be
                # prescribed

                else:

                    if isinstance(load_info[0], list):

                        # Iterates through the prescribed DOFs

                        for dof_number in load_info[0]:

                            # Adds this particular boundary condition 
                            # to the lot

                            boundary_conditions.append(DirichletBC(
                            field_functionSpace.sub(dof_number), 
                            load_info[1], boundary_meshFunction, 
                            physical_group))

                    else:

                        # Adds this particular boundary condition to 
                        # the lot

                        boundary_conditions.append(DirichletBC(
                        field_functionSpace.sub(load_info[0]), 
                        load_info[1], boundary_meshFunction, 
                        physical_group))

            # Otherwise, iterates through the fields

            else:

                for field in range(n_fields):

                    # Verifies if this field is to be constrained or 
                    # not

                    if field==load_info[0]:

                        # Checks whether a set of degrees of freedom
                        # of the field is to be prescribed

                        if len(load_info)>2:

                            # Verifies if there is a list of DOFs

                            if isinstance(load_info[1], list):

                                # Iterates through the prescribed
                                # DOFs

                                for dof_number in load_info[1]:

                                    # Adds this particular boundary
                                    # condition to the lot

                                    boundary_conditions.append(DirichletBC(
                                    field_functionSpace.sub(field).sub(
                                    dof_number), load_info[2], 
                                    boundary_meshFunction, 
                                    physical_group))

                            else:

                                # Adds this particular boundary con-
                                # dition to the lot

                                boundary_conditions.append(DirichletBC(
                                field_functionSpace.sub(field).sub(
                                load_info[1]), load_info[2], 
                                boundary_meshFunction, 
                                physical_group))

                        else:

                            # Adds this particular boundary condi-
                            # tion to the lot

                            boundary_conditions.append(DirichletBC(
                            field_functionSpace.sub(field), 
                            load_info[1], boundary_meshFunction, 
                            physical_group))

    # Returns the boundary conditions list

    if verbose:

        print("Finishes creating fixed support boundary conditions.\n")

    return boundary_conditions

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
boundary_physicalGroups=0, sub_fieldsToApplyBC=[], boundary_conditions=[
], boundary_physGroupsNamesToTags=dict(), verbose=False):

    # If the physical groups variable is null, returns the empty list of
    # boundary conditions

    if boundary_physicalGroups==0:

        if verbose:

            print("Creates no fixed support boundary condition.\n")

        return boundary_conditions
    
    # Tests if the element is mixed and gets the number of fields

    n_fields = 1

    element_type = field_functionSpace.ufl_element().family()

    if element_type=='Mixed':

        # Gets the number of fields

        n_fields = field_functionSpace.ufl_element().num_sub_elements()

    # Verifies whether the boundary physical groups is a list or not

    if isinstance(boundary_physicalGroups, list):

        # Iterates through the regions

        for physical_group in boundary_physicalGroups:

            # Verifies if the physical group is a string

            physical_group = verify_stringPhysicalGroup(physical_group, 
            boundary_physGroupsNamesToTags)

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

        # Verifies if the physical group is a string

        boundary_physicalGroups = verify_stringPhysicalGroup(
        boundary_physicalGroups, boundary_physGroupsNamesToTags)

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

    if verbose:

        print("Finishes creating fixed support boundary conditions.\n")

    return boundary_conditions

# Defines a function to apply a simple support. If sub_fields is an empty
# list, all DOFs are constrained

def simple_supportDirichletBC(field_functionSpace, boundary_meshFunction, 
boundary_physicalGroups, sub_fieldsToApplyBC=[], boundary_conditions=[],
boundary_physGroupsNamesToTags=dict(), verbose=False):
    
    # Verifies if the boundary physical groups is a dictionary

    if not isinstance(boundary_physicalGroups, dict):

        raise ValueError("The boundary_physicalGroups variable in the "+
        "simple_supportDirichletBC method must be a dictionary, where "+
        "the keys are the physical groups of the boundary regions and "+
        "values are the list of DOFs to be constrained.")

    # If the physical groups variable is null, returns the empty list of
    # boundary conditions

    if len(list(boundary_physicalGroups.keys()))==0:

        if verbose:

            print("Creates no simple support boundary condition.\n")

        return boundary_conditions
    
    # Tests if the element is mixed and gets the number of fields

    n_fields = 1

    element_type = field_functionSpace.ufl_element().family()

    if element_type=='Mixed':

        # Gets the number of fields

        n_fields = field_functionSpace.ufl_element().num_sub_elements()

    # Iterates through the regions

    for physical_group, list_constrainedDOFs in boundary_physicalGroups.items():

        # Verifies if the physical group is a string

        physical_group = verify_stringPhysicalGroup(physical_group, 
        boundary_physGroupsNamesToTags)

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

    # Returns the boundary conditions list

    if verbose:

        print("Finishes creating fixed support boundary conditions.\n")

    return boundary_conditions

########################################################################
#                              Utilities                               #
########################################################################

# Defines a function to verify if the physical group is a string and, 
# then, to convert it into a number tag using the dictionary of physical
# groups' names to tags

def verify_stringPhysicalGroup(physical_group, 
boundary_physGroupsNamesToTags):
    
    if isinstance(physical_group, str):

        # Converts this physical group to the corresponding number tag

        try:

            physical_group = boundary_physGroupsNamesToTags[
            physical_group]

        except:

            raise KeyError("The physical group '"+physical_group+"' is"+
            " not in the dictionary of physical groups names' to tags."+
            " This dictionary has the following keys and values: "+str(
            boundary_physGroupsNamesToTags))
        
    return physical_group