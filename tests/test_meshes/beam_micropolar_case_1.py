# Routine to store the generation of the beam case of the micropolar 
# multiscale analysis paper

import source.solids.cuboid_cylinders as cylinders

import source.tool_box.meshing_tools as tools

import source.tool_box.region_finder as region_finder

def case_1(RVE_width, RVE_length, fiber_radius, n_RVEsX, n_RVEsY, 
n_RVEsZ, RVE_localizationX, RVE_localizationY, RVE_localizationZ, 
mesh_fileName="micropolar_beam_with_fibers", file_directory=
"tests//micropolar_meshes//results"):

    ####################################################################
    #                     RVE geometric properties                     #
    ####################################################################

    # Defines the longitudinal direction of the fiber

    fiber_direction = [1.0, 0.0, 0.0]

    # Defines the normal vectors and the biases of the planes that cons-
    # traint the RVE

    RVE_planesNormals = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 
    -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

    RVE_planesBiases = [-(RVE_localizationX-1)*RVE_length, (-(
    RVE_localizationY-1)*RVE_width), -(RVE_localizationZ-1)*RVE_width, (
    RVE_localizationX*RVE_length), RVE_localizationY*RVE_width,
    RVE_localizationZ*RVE_width]

    ####################################################################
    #                          Mesh parameters                         #
    ####################################################################

    # Defines the mesh discretization and biases

    transfinite_directions = [6, 6, 3, 4, 3]

    bias_directions = dict()

    bias_directions["x"] = 1.0

    bias_directions["y"] = 1.0

    bias_directions["z"] = 1.0

    bias_directions["cylinder radial"] = 1.5

    bias_directions["box radial"] = 1.5

    ####################################################################
    #                         Physical groups                          #
    ####################################################################

    # Sets the finder for the boundary surfaces

    # XY plane at z = 0

    def back_expression(x, y, z):

        return z

    # XZ plane at y = 0

    def bottom_expression(x, y, z):

        return y

    # XZ plane at y = n_RVEsY*RVE_width

    def top_expression(x, y, z):

        return y-(n_RVEsY*RVE_width)

    # Sets a list of expressions to find the surfaces at the boundaries

    surface_regionsExpressions = [back_expression, bottom_expression, 
    top_expression]

    # Sets the names of the surface regions

    surface_regionsNames = ['back', 'lower', 'upper']

    # Sets and expression to find any fiber testing different fibers

    def find_fiber(x, y, z):

        # Iterates through the arrangement of fiber

        for i in range(n_RVEsY):

            for j in range(n_RVEsZ):

                # Verifies if this point is in this particular fiber and
                # not inside and RVE
                
                if (region_finder.cylindrical_enclosure(x, y, z, 
                fiber_radius, RVE_length*n_RVEsX, fiber_direction, [[
                -fiber_direction[0], -fiber_direction[1], 
                -fiber_direction[2]], fiber_direction], [0.0, ((i+0.5)*
                RVE_width), ((j+0.5)*RVE_width)], flag_inside=True, 
                tolerance=1E-4)):
                    
                    return True
                
        # If this point was not found in a fiber, returns False

        return False

    # Sets and expression to find the matrix

    def find_matrix(x, y, z):

        # Iterates through the arrangement of fiber

        for i in range(n_RVEsY):

            for j in range(n_RVEsZ):

                # Verifies if this point is in this particular fiber
                
                if (not region_finder.cylindrical_enclosure(x, y, z, 
                fiber_radius, RVE_length*n_RVEsX, fiber_direction, [[
                -fiber_direction[0], -fiber_direction[1], 
                -fiber_direction[2]], fiber_direction], [0.0, ((i+0.5)*
                RVE_width), ((j+0.5)*RVE_width)], flag_inside=False, 
                tolerance=1E-4)):
                    
                    return False
                
        # If this point was not found in a fiber, returns True, because 
        # it must be in the matrix

        return True

    # Sets and expression to find the RVE fiber testing different fibers

    def find_RVEfiber(x, y, z):

        if (region_finder.cylindrical_enclosure(x, y, z, fiber_radius, 
        RVE_length*n_RVEsX, fiber_direction, [[-fiber_direction[0], 
        -fiber_direction[1], -fiber_direction[2]], fiber_direction], [
        0.0, (RVE_localizationY-0.5)*RVE_width, ((RVE_localizationZ-0.5)
        *RVE_width)], flag_inside=True, tolerance=1E-4) and (
        region_finder.plane_enclosure(x, y, z, RVE_planesNormals, 
        RVE_planesBiases))):
                    
            return True

        return False

    # Sets and expression to find the RVE matrix

    def find_RVEmatrix(x, y, z):

        if (region_finder.cylindrical_enclosure(x, y, z, fiber_radius, 
        RVE_length*n_RVEsX, fiber_direction, [[-fiber_direction[0], 
        -fiber_direction[1], -fiber_direction[2]], fiber_direction], [
        0.0, (RVE_localizationY-0.5)*RVE_width, ((RVE_localizationZ-0.5)
        *RVE_width)], flag_inside=False, tolerance=1E-4) and (
        region_finder.plane_enclosure(x, y, z, RVE_planesNormals, 
        RVE_planesBiases))):
                    
            return True
        
        return False

    # Sets the list of volume expressions

    volume_regionsExpressions = [find_RVEfiber, find_RVEmatrix, 
    find_fiber, find_matrix]

    # Sets the list of volume physical groups' names

    volume_regionsNames = ['RVE fiber', 'RVE matrix', 'Fiber', 'Matrix']

    ####################################################################
    #                       Geometry generation                        #
    ####################################################################

    # Initializes the gmsh 

    geometric_data = tools.gmsh_initialization(
    surface_regionsExpressions=surface_regionsExpressions, 
    volume_regionsExpressions=volume_regionsExpressions, 
    surface_regionsNames=surface_regionsNames, volume_regionsNames=
    volume_regionsNames)

    # Creates the RVEs iterating though the directions of reptibility of 
    # RVEs

    for i in range(n_RVEsX):

        for j in range(n_RVEsY):

            for k in range(n_RVEsZ):
                
                # Creates the base point where the RVE is centered, the 
                # centroid of the lower facet

                base_point = [i*RVE_length, (j+0.5)*RVE_width, ((k+0.5)*
                RVE_width)]

                geometric_data = cylinders.cylinder_inBox(fiber_radius, 
                fiber_radius, RVE_length, RVE_width, RVE_width, 
                fiber_direction, fiber_direction, fiber_direction, 
                base_point, transfinite_directions=
                transfinite_directions, bias_directions=bias_directions, 
                geometric_data=geometric_data)

    ####################################################################
    #                 Mesh generation and file creation                #
    ####################################################################

    # Finalizes the mesh and saves the result

    tools.gmsh_finalize(geometric_data=geometric_data, file_directory=
    file_directory, file_name=mesh_fileName)