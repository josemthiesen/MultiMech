# Routine to store methods to neatly apply boundary conditions

from dolfin import *

import numpy as np

import source.tool_box.programming_tools as programming_tools

import source.tool_box.mesh_handling_tools as mesh_tools

import source.tool_box.numerical_tools as numerical_tools

########################################################################
#              Heterogeneous Dirichlet boundary conditions             #
########################################################################

# Defines a function to apply prescribed boundary conditions to a pro-
# blem using an Expression as load prescriber

@programming_tools.optional_argumentsInitializer({'boundary_conditions':
lambda: []})

def PrescribedDirichletBC(prescribed_conditionsDict, 
field_functionSpace, mesh_dataClass, fields_namesDict, 
boundary_conditions=None, dirichlet_loads=None, t_initial=0.0):
    
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

        # Verifies if load info is a list

        if not isinstance(load_info, list):

            raise TypeError("The load_info for the prescribed Dirichle"+
            "t boundary condition is not a list, but: "+str(load_info)+
            "\nIt must be a list")
        
        # Verifies if the load itself is a string, which denotes that it
        # comes from a loading curve

        elif isinstance(load_info[-1], str):

            # Initializes the initial time constant

            time_constant = Constant(t_initial)

            # Gets the loading curve

            load_info[-1] = numerical_tools.generate_loadingParametricCurves(
            load_info[-1])(time_constant)

            # Adds the time constant to the dirichlet_loads

            dirichlet_loads.append(time_constant)

        # Verifies whether the boundary physical groups is a list or not

        if isinstance(physical_groups, tuple):

            # Iterates through the regions

            for physical_group in physical_groups:

                # Verifies if the physical group is a string

                physical_group = verify_stringPhysicalGroup(physical_group, 
                mesh_dataClass.boundary_physicalGroupsNameToTag)

                # Verifies if there is only one field

                if n_fields==1:

                    # Verifies if no DOFs are especified

                    if len(load_info)==1:

                        # Adds this particular boundary condition to the 
                        # lot

                        boundary_conditions.append(DirichletBC(
                        field_functionSpace, load_info[0], 
                        mesh_dataClass.boundary_meshFunction, 
                        physical_group))

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
                                load_info[1], 
                                mesh_dataClass.boundary_meshFunction, 
                                physical_group))

                        else:

                            # Adds this particular boundary condition to 
                            # the lot

                            boundary_conditions.append(DirichletBC(
                            field_functionSpace.sub(load_info[0]), 
                            load_info[1], 
                            mesh_dataClass.boundary_meshFunction, 
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
                                        mesh_dataClass.boundary_meshFunction, 
                                        physical_group))

                                else:

                                    # Adds this particular boundary con-
                                    # dition to the lot

                                    boundary_conditions.append(DirichletBC(
                                    field_functionSpace.sub(field).sub(
                                    load_info[1]), load_info[2], 
                                    mesh_dataClass.boundary_meshFunction, 
                                    physical_group))

                            else:

                                # Adds this particular boundary condi-
                                # tion to the lot

                                boundary_conditions.append(DirichletBC(
                                field_functionSpace.sub(field), 
                                load_info[1], 
                                mesh_dataClass.boundary_meshFunction, 
                                physical_group))

        # Otherwise, if there is only one physical group to apply boundary 
        # conditions

        else:

            # Verifies if the physical group is a string

            physical_group = verify_stringPhysicalGroup(physical_groups, 
            mesh_dataClass.boundary_physicalGroupsNameToTag)

            # Verifies if there is only one field

            if n_fields==1:

                # Verifies if no DOFs are especified

                if len(load_info)==1:

                    # Adds this particular boundary condition to the lot

                    boundary_conditions.append(DirichletBC(
                    field_functionSpace, load_info[0], 
                    mesh_dataClass.boundary_meshFunction, physical_group
                    ))

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
                            load_info[1], 
                            mesh_dataClass.boundary_meshFunction, 
                            physical_group))

                    else:

                        # Adds this particular boundary condition to 
                        # the lot

                        boundary_conditions.append(DirichletBC(
                        field_functionSpace.sub(load_info[0]), 
                        load_info[1], 
                        mesh_dataClass.boundary_meshFunction, 
                        physical_group))

            # Otherwise, iterates through the fields

            else:

                # Verifies if the load info has string values

                load_info[0] = convert_fieldsNamesToIndices(load_info[0], 
                fields_namesDict, dont_convertToList=True)

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
                                    mesh_dataClass.boundary_meshFunction, 
                                    physical_group))

                            else:

                                # Adds this particular boundary con-
                                # dition to the lot

                                boundary_conditions.append(DirichletBC(
                                field_functionSpace.sub(field).sub(
                                load_info[1]), load_info[2], 
                                mesh_dataClass.boundary_meshFunction, 
                                physical_group))

                        else:

                            # Adds this particular boundary condi-
                            # tion to the lot

                            boundary_conditions.append(DirichletBC(
                            field_functionSpace.sub(field), 
                            load_info[1], 
                            mesh_dataClass.boundary_meshFunction, 
                            physical_group))

    # Returns the boundary conditions list

    if mesh_dataClass.verbose:

        print("Finishes creating prescribed Dirichlet boundary conditi"+
        "ons.\n")

    if dirichlet_loads is None:

        return boundary_conditions
    
    else:

        return boundary_conditions, dirichlet_loads

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

@programming_tools.optional_argumentsInitializer({'boundary_conditions':
lambda: [], 'sub_fieldsToApplyBC': lambda: []})

def FixedSupportDirichletBC(field_functionSpace, mesh_dataClass,
fields_namesDict, boundary_physicalGroups=0, sub_fieldsToApplyBC=None, 
boundary_conditions=None):

    # If the physical groups variable is null, returns the empty list of
    # boundary conditions

    if boundary_physicalGroups==0:

        if mesh_dataClass.verbose:

            print("Creates no fixed support boundary condition.\n")

        return boundary_conditions
    
    # If subfields are required, use the dictionary of fields to convert
    # the names to the numbers

    if not (sub_fieldsToApplyBC is None):

        sub_fieldsToApplyBC = convert_fieldsNamesToIndices(
        sub_fieldsToApplyBC, fields_namesDict)
    
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
            mesh_dataClass.boundary_physicalGroupsNameToTag)

            # Verifies if there is only one field

            if n_fields==1:

                # Adds this particular boundary condition to the lot

                boundary_conditions.append(DirichletBC(
                field_functionSpace, Constant((0.0, 0.0, 0.0)), 
                mesh_dataClass.boundary_meshFunction, physical_group))

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
                        0.0)), mesh_dataClass.boundary_meshFunction, 
                        physical_group))

    # Otherwise, if there is only one physical group to apply boundary 
    # conditions

    else:

        # Verifies if the physical group is a string

        boundary_physicalGroups = verify_stringPhysicalGroup(
        boundary_physicalGroups, 
        mesh_dataClass.boundary_physicalGroupsNameToTag)

        # Verifies if there is only one field

        if n_fields==1:

            # Adds this particular boundary condition to the lot

            boundary_conditions.append(DirichletBC(field_functionSpace,
            Constant((0.0, 0.0, 0.0)), 
            mesh_dataClass.boundary_meshFunction, boundary_physicalGroups
            ))

        # Otherwise, iterates through the fields

        else:

            for field in range(n_fields):

                # Verifies if this field is to be constrained or not

                if (field in sub_fieldsToApplyBC or (
                sub_fieldsToApplyBC==[])):

                    # Adds this particular boundary condition to the lot

                    boundary_conditions.append(DirichletBC(
                    field_functionSpace.sub(field), Constant((0.0, 0.0, 
                    0.0)), mesh_dataClass.boundary_meshFunction, 
                    boundary_physicalGroups))

    # Returns the boundary conditions list

    if mesh_dataClass.verbose:

        print("Finishes creating fixed support boundary conditions.\n")

    return boundary_conditions

# Defines a function to apply a simple support. If sub_fields is an empty
# list, all DOFs are constrained

@programming_tools.optional_argumentsInitializer({'boundary_conditions':
lambda: [], 'sub_fieldsToApplyBC': lambda: []})

def SimpleSupportDirichletBC(field_functionSpace, mesh_dataClass, 
boundary_physicalGroupsDict, fields_namesDict, sub_fieldsToApplyBC=None, 
boundary_conditions=None, boundary_physicalGroups=None):
    
    # Verifies if the boundary physical groups is a dictionary or a list
    # if the physical group is given

    if isinstance(boundary_physicalGroupsDict, list):

        if boundary_physicalGroups is None:

            raise ValueError("To apply SimpleSupportDirichletBC, the b"+
            "oundary_physicalGroups variables can be a list if and onl"+
            "y if the boundary_physicalGroups is provided, but it's No"+
            "ne")
        
        # Transforms in a dictionary
        
        boundary_physicalGroupsDict[boundary_physicalGroups
        ] = boundary_physicalGroupsDict

    elif not isinstance(boundary_physicalGroupsDict, dict):

        raise ValueError("The boundary_physicalGroupsDict variable in the "+
        "simple_supportDirichletBC method must be a dictionary, where "+
        "the keys are the physical groups of the boundary regions and "+
        "values are the list of DOFs to be constrained.")

    # If the physical groups variable is null, returns the empty list of
    # boundary conditions

    if len(list(boundary_physicalGroupsDict.keys()))==0:

        if mesh_dataClass.verbose:

            print("Creates no simple support boundary condition.\n")

        return boundary_conditions
    
    # If subfields are required, use the dictionary of fields to convert
    # the names to the numbers

    if not (sub_fieldsToApplyBC is None):

        sub_fieldsToApplyBC = convert_fieldsNamesToIndices(
        sub_fieldsToApplyBC, fields_namesDict)
    
    # Tests if the element is mixed and gets the number of fields

    n_fields = 1

    element_type = field_functionSpace.ufl_element().family()

    if element_type=='Mixed':

        # Gets the number of fields

        n_fields = field_functionSpace.ufl_element().num_sub_elements()

    # Iterates through the regions

    for physical_group, list_constrainedDOFs in boundary_physicalGroupsDict.items():

        # Verifies if the physical group is a string

        physical_group = verify_stringPhysicalGroup(physical_group, 
        mesh_dataClass.boundary_physicalGroupsDictNameToTag)

        # Verifies if there is only one field

        if n_fields==1:

            # Iterates through the DOFs to be constrained

            for i in range(field_functionSpace.num_sub_spaces()):
                
                if i in list_constrainedDOFs:

                    # Adds this particular boundary condition to the 
                    # lot

                    boundary_conditions.append(DirichletBC(
                    field_functionSpace.sub(i), Constant(0.0), 
                    mesh_dataClass.boundary_meshFunction, physical_group
                    ))

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
                            Constant(0.0), 
                            mesh_dataClass.boundary_meshFunction, 
                            physical_group))

    # Returns the boundary conditions list

    if mesh_dataClass.verbose:

        print("Finishes creating simply supported boundary conditions."+
        "\n")

    return boundary_conditions

########################################################################
#                          Subdomain classes                           #
########################################################################

# Defines a function to generate a class to apply boundary conditions to 
# a node

def generate_nodeSubdomain(point_coordinates, mesh_dataClass, tolerance=
1E-5):

    # Gets the coordinates of the node closest to the point required

    _, node_coordinates = mesh_tools.find_nodeClosestToPoint(
    mesh_dataClass, point_coordinates, None, None)

    # Defines the class

    class FixedNode(SubDomain):

        def __init__(self, node_coordinates, tolerance):

            super().__init__()

            self.node_coordinates = node_coordinates

            self.tolerance = tolerance

        def inside(self, x, on_boundary):

            # Ignores the on_boudnary flag because a node inside the do-
            # main can be fixed

            return (near(x[0], self.node_coordinates[0], self.tolerance
            ) and near(x[1], self.node_coordinates[1], self.tolerance)
            and near(x[2], self.node_coordinates[2], self.tolerance))

    # Instantiates the class

    fixed_node = FixedNode(node_coordinates, tolerance)

    return fixed_node

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

# Defines a function to convert a list or a value of fields to a list of
# their corresponding indices

def convert_fieldsNamesToIndices(sub_fieldsToApplyBC, fields_namesDict,
dont_convertToList=False):

    if isinstance(sub_fieldsToApplyBC, list):

        for i in range(len(sub_fieldsToApplyBC)):

            if isinstance(sub_fieldsToApplyBC[i], str):

                if sub_fieldsToApplyBC[i] in fields_namesDict:

                    sub_fieldsToApplyBC[i] = fields_namesDict[
                    sub_fieldsToApplyBC[i]]

                else:

                    raise KeyError("The field '"+str(sub_fieldsToApplyBC[
                    i])+"' was not found in the dictionary of fields' "+
                    "names. Hence, it is not possible to make a fixed "+
                    "support for this field. Check out the list of ava"+
                    "ilable fields names: "+str(fields_namesDict.keys()))
                
            elif not isinstance(sub_fieldsToApplyBC[i], int):

                # If the field is not an integer, the index of the field,
                # throws an error

                raise ValueError("Only integer or strings can be used "+
                "to recover fields to make fixed supports. The value "+
                str(sub_fieldsToApplyBC[i])+" was given")
            
    elif isinstance(sub_fieldsToApplyBC, str):
        
        if sub_fieldsToApplyBC in fields_namesDict:

            if dont_convertToList:

                # Does not transform to a list

                sub_fieldsToApplyBC = fields_namesDict[
                sub_fieldsToApplyBC]

            else:

                # Transforms to a list also

                sub_fieldsToApplyBC = [fields_namesDict[
                sub_fieldsToApplyBC]]

        else:

            raise KeyError("The field '"+str(sub_fieldsToApplyBC)+"' w"+
            "as not found in the dictionary of fields' names. Hence, i"+
            "t is not possible to make a fixed support for this field."+
            " Check out the list of available fields names: "+str(
            fields_namesDict.keys()))
        
    elif not isinstance(sub_fieldsToApplyBC, int):

        # If the field is not an integer, the index of the field, throws
        # an error

        raise ValueError("Only integer or strings can be used to recov"+
        "er fields to make fixed supports. The value "+str(
        sub_fieldsToApplyBC[i])+" was given")
    
    return sub_fieldsToApplyBC