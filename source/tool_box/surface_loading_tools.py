# Routine to store tools to create loads

from dolfin import *

import ufl_legacy as ufl

import numpy as np

import source.tool_box.mesh_handling_tools as mesh_tools

########################################################################
#                          Referential forces                          #
########################################################################

# Defines a function to construct a uniform traction on a surface in the
# referential configuration

def UniformReferentialTraction(amplitude_tractionX, amplitude_tractionY,
amplitude_tractionZ, t=0.0, t_final=1.0, parametric_load_curve=lambda x: 
x):

    # Creates the time constants

    time_constant = Constant(t)

    maximum_time = Constant(t_final)

    # Evaluates the traction vector at the referential configuration

    T = as_vector([parametric_load_curve(time_constant/maximum_time)*
    amplitude_tractionX, parametric_load_curve(time_constant/
    maximum_time)*amplitude_tractionY, parametric_load_curve(
    time_constant/maximum_time)*amplitude_tractionZ])

    # Returns the traction and the time constant

    return T, time_constant

########################################################################
#                           Follower forces                            #
########################################################################

# Defines a function to construct a normal uniform traction in a surfa-
# ce, this means the traction vector is always perpendicular to the sur-
# face. The user can set how the amplitude of the function is interpola-
# ted between the initial and final time points; the default value is a
# linear curve

def NormalUniformFollowerTraction(field, mesh_dataClass,
amplitude_traction, t=0.0, t_final=1.0, parametric_load_curve=lambda x: 
x):

    # Creates the time constants

    time_constant = Constant(t)

    maximum_time = Constant(t_final)

    # Evaluates the deformation gradient

    I = Identity(3)

    F = grad(field)+I

    J = det(F)

    # Evaluates the traction vector at the referential configuration

    T = (parametric_load_curve(time_constant/maximum_time)*
    amplitude_traction*J*((inv(F)).T)*mesh_dataClass.n)

    # Returns the traction and the time constant

    return T, time_constant

# Defines a function to construct a torsion on a surface. The center 
# point, where the axis of torsion intersect the surface, is optional,
# if not given, the centroid of the surface will be used automatically.
# The torsion is applied as a shearing traction inside a circular re-
# gion, whose center is the center point; the radius can be given or e-
# valuated by the minimum distance to the edge if not supplied

def NormalFollowerTorsion(field, mesh_dataClass, amplitude_torsion, 
physical_group, center_point=None, influence_radius=None, t=0.0, t_final
=1.0, parametric_load_curve=lambda x: x):
    
    # Verifies if the physical group is a string

    if isinstance(physical_group, str):

        # Tests if it is in the dictionary of physical groups of the 
        # boundary

        if physical_group in mesh_dataClass.boundary_physicalGroupsNameToTag:

            # Converts it

            physical_group = mesh_dataClass.boundary_physicalGroupsNameToTag[
            physical_group]

        else:

            raise KeyError("The physical group '"+str(physical_group)+
            "' is not in the dictionary of boundary physical groups. T"+
            "hus, cannot be used to find the nodes in the boundary of "+
            "surface. Check out the available options of physical grou"+
            "ps in the boundary: "+str(
            mesh_dataClass.boundary_physicalGroupsNameToTag.keys()))
        
    # Does not accept any other formats than integer

    elif not isinstance(physical_group, int):

        raise TypeError("The physical group "+str(physical_group)+" is"+
        " not an integer, thus cannot be used to find the nodes in the"+
        " boundary of a surface")

    # Verifies if the center point is None, then, proceeds to evaluate
    # the centroid of the region

    if center_point is None:

        # Gets the position vector from the mesh

        position_vector = mesh_dataClass.x

        # Evaluates the area of this physical group

        area_inverse = (1.0/float(assemble(1.0*mesh_dataClass.ds(
        physical_group))))

        # Evaluates the centroid coordinates

        centroid_x = (area_inverse*float(assemble(position_vector[0]*ds(
        physical_group))))

        centroid_y = (area_inverse*float(assemble(position_vector[1]*ds(
        physical_group))))

        centroid_z = (area_inverse*float(assemble(position_vector[2]*ds(
        physical_group))))

        center_point = [centroid_x, centroid_y, centroid_z]

    # Otherwise, verifies if it is a list

    elif not (isinstance(center_point, list)):

        raise TypeError("The center point of application of the follow"+
        "er torsion must be a list. The following was given though: "+
        str(center_point))

    # Gets the node closest to the center point

    node_number, center_point = mesh_tools.find_nodeClosestToPoint(
    mesh_dataClass, center_point, None, None)

    # If the radius of influence was not given, evaluate it as the mini-
    # mum distance from the center point to the edge of the surface

    if influence_radius is None:

        # Gets the nodes on the boundary of the surface

        boundary_nodesSet = mesh_tools.find_nodesOnSurfaceBoundary(
        mesh_dataClass, physical_group)

        # Gets the node number and its coordinates where is closest to
        # the center point

        node_number, node_coordinates = mesh_tools.find_nodeClosestToPoint(
        mesh_dataClass, center_point, None, None, set_ofNodes=
        boundary_nodesSet)

        # Gets the radius 

        influence_radius = np.sqrt(((node_coordinates[0]-center_point[0]
        )**2)+((node_coordinates[1]-center_point[1])**2)+((
        node_coordinates[2]-center_point[2])**2))

    # Otherwise, checks if it is a float

    elif not (isinstance(influence_radius, int) or isinstance(
    influence_radius, float)):
        
        raise TypeError("The influence radius of the NormalFollowerTor"+
        "sion must be a float number, whereas the given value is "+str(
        influence_radius))
    
    # Updates the amplitude torsion to be the linear coefficient of the
    # radially linear distribution of shearing traction. This ensures 
    # that the amplitude torsion is indeed the moment applied to the 
    # surface

    #amplitude_torsion = ((2.0*amplitude_torsion)/(np.pi*(
    #influence_radius**4)))

    amplitude_torsion = Constant(amplitude_torsion)

    # Creates the time constants

    time_constant = Constant(t)

    maximum_time = Constant(t_final)
    
    # Evaluates the deformation gradient

    I = Identity(3)

    F = grad(field)+I

    J = det(F)

    # Gets the normal vector in the deformed configuration (but it does
    # not have to be unitary because, later, the whole vector will be 
    # rescaled)

    normal_vector = (inv(F).T)*mesh_dataClass.n

    # Sets the centroid vector as constant

    centroid_point = Constant(center_point)

    # Sets the relative position vector

    relative_position = mesh_dataClass.x-centroid_point

    # Gets the relative displacement (subtracting the displacement at 
    # the center point)

    centroid_displacement = field(Point(center_point))

    relative_displacement = as_vector([field[0]-centroid_displacement[0],
    field[1]-centroid_displacement[1], field[2]-centroid_displacement[2]
    ])

    # Gets the traction direction as the cross product of the relative
    # position translated by the displacement with the normal vector

    traction_direction = cross(normal_vector, relative_position+
    relative_displacement)

    # Gets the norm of this vector. Multiplies by the magnitude of the 
    # relative position vector squared and divided by the influence ra-
    # dius squared to get a radially increasing traction profile

    norm_traction = (dot(relative_position, relative_position)/
    ufl.conditional(gt(sqrt(dot(traction_direction, traction_direction)
    ), 1E-5), sqrt(dot(traction_direction, traction_direction)), 1.0))

    # Creates the traction vector. Multiplies by a correction factor to
    # account for the pull back to the referential configuration

    T_unscaled = (norm_traction*traction_direction*J*sqrt(dot(normal_vector, 
    normal_vector)))

    # Gets the average normal vector of the facet

    average_normalVector = [(area_inverse*normal_vector[0]*
    mesh_dataClass.ds(physical_group)), (area_inverse*normal_vector[1]*
    mesh_dataClass.ds(physical_group)), (area_inverse*normal_vector[2]*
    mesh_dataClass.ds(physical_group))]

    # Evaluates the resulting moment on the facet

    d_moment = cross(relative_position, T_unscaled)

    moment = (sqrt(dot(d_moment, d_moment))*mesh_dataClass.ds(
    physical_group))

    # Scales the traction vector to give the actual moment

    T = (parametric_load_curve(time_constant/maximum_time)*(
    amplitude_torsion/moment)*T_unscaled)

    return T, time_constant