# Routine to store user expressions for multiscale boundary conditions

from dolfin import *

import numpy as np

########################################################################
#                   Boundary conditions' expressions                   #
########################################################################

# Defines a class for an expression for a scalar linear field with res-
# pect to space coordinates

class LinearScalarFieldExpression(UserExpression):

    def __init__(self, mean_field, gradient_field, x_centroid, 
    y_centroid, z_centroid):

        super().__init__()

        # Initializes the values

        self.mean_field = mean_field

        self.gradient_field = gradient_field

        self.x_centroid = x_centroid

        self.y_centroid = y_centroid

        self.z_centroid = z_centroid

    # Defines a function to evaluate the expression

    def eval(self, value, x):

        # Gets the position vector shifting the coordinates vector by 
        # the centroid

        y = [x[0]-self.x_centroid, x[1]-self.y_centroid, (x[2]-
        self.z_centroid)]

        # Checks if the mean field is a list or else, if it is a Cons-
        # tant. If it is a Constant, gets the numerical values

        if isinstance(self.mean_field, list):

            mean_field = self.mean_field

        else:

            mean_field = [a for a in self.mean_field.values()]

        # Checks if the gradient of field is a list or else, if it is a
        # Constant. If it is a Constant, gets the numerical values

        if isinstance(self.gradient_field, list):

            gradient_field = self.gradient_field

        else:

            gradient_field = self.gradient_field.values()

        # Updates the values of the expression to the linear profile 
        # that is required

        value = (mean_field[0]+(gradient_field[0]*y[0])+(
        gradient_field[1]*y[1])+(gradient_field[2]*y[2]))

    def value_shape(self):

        return ()

# Defines a class for an expression for a vector linear field with res-
# pect to space coordinates

class LinearVectorFieldExpression(UserExpression):

    def __init__(self, mean_field, gradient_field, x_centroid, 
    y_centroid, z_centroid):

        super().__init__()

        # Initializes the values

        self.mean_field = mean_field

        self.gradient_field = gradient_field

        self.x_centroid = x_centroid
        
        self.y_centroid = y_centroid
        
        self.z_centroid = z_centroid

    # Defines a function to evaluate the expression

    def eval(self, value, x):

        # Gets the position vector shifting the coordinates vector by 
        # the centroid

        y = [x[0]-self.x_centroid, x[1]-self.y_centroid, (x[2]-
        self.z_centroid)]

        # Checks if the mean field is a list or else, if it is a Cons-
        # tant. If it is a Constant, gets the numerical values

        if isinstance(self.mean_field, list):

            mean_field = self.mean_field

        else:

            mean_field = [a for a in self.mean_field.values()]

        # Checks if the gradient of field is a list or else, if it is a
        # Constant. If it is a Constant, gets the numerical values

        if isinstance(self.gradient_field, list):

            gradient_field = self.gradient_field

        else:

            A = self.gradient_field.values()

            gradient_field = [[A[0], A[1], A[2]], [A[3], A[4], A[5]], [
            A[6], A[7], A[8]]]

        # Updates the values of the expression to the linear profile 
        # that is required

        value[0] = (mean_field[0]+(gradient_field[0][0]*y[0])+
        (gradient_field[0][1]*y[1])+(gradient_field[0][2]*y[2]
        ))

        value[1] = (mean_field[1]+(gradient_field[1][0]*y[0])+
        (gradient_field[1][1]*y[1])+(gradient_field[1][2]*y[2]
        ))

        value[2] = (mean_field[2]+(gradient_field[2][0]*y[0])+
        (gradient_field[2][1]*y[1])+(gradient_field[2][2]*y[2]
        ))

    def value_shape(self):

        return (3,)

########################################################################
#                  Boundary conditions' domain classes                 #
########################################################################

# Defines a class for a function space to have a domain constrained as
# periodic. The domain must be a cube (the dimensions don't have to be
# the same)

class PeriodicBoundary(Subdomain):

    def __init__(self, mean_field, gradient_field, x_centroid, 
    y_centroid, z_centroid, length_x, length_y, length_z):

        super().__init__()

        # Initializes the values

        self.mean_field = mean_field

        self.gradient_field = gradient_field

        self.x_centroid = x_centroid
        
        self.y_centroid = y_centroid
        
        self.z_centroid = z_centroid

        # Gets the cube dimensions

        self.semi_lengthX = 0.5*length_x

        self.semi_lengthY = 0.5*length_y

        self.semi_lengthZ = 0.5*length_z

    # Redefines the inside method to give all True for points on the 
    # master facets, which are conveined here as x=-0.5*length_x, y=-0.5
    # *length_y, and z=-0.5*length_z planes. The slave facets are falsi-
    # fied because they are dependent on the master ones

    def inside(self, x, on_boundary):

        flag_xFacet = near(self.x_centroid-x[0], self.semi_lengthX)

        flag_yFacet = near(self.y_centroid-x[1], self.semi_lengthY)

        flag_zFacet = near(self.z_centroid-x[2], self.semi_lengthZ)

        return bool ((flag_xFacet or flag_yFacet or flag_zFacet) and
        on_boundary)

    # Defines a method to map the field from one facet to the opposite

    def map(self, x, y):

        # If the point x is on the vertex of the three facets, maps it 
        # to the opposite point

        if (near(x[0], self.x_centroid-self.semi_lengthX) and near(x[1], 
        self.y_centroid-self.semi_lengthY) and near(x[2], 
        self.z_centroid-self.semi_lengthZ)):

            y[0] = x[0]+(2*self.semi_lengthX)

            y[1] = x[1]+(2*self.semi_lengthY)

            y[2] = x[2]+(2*self.semi_lengthZ)

        # If the point x is in the edge along the X axis

        if (near(x[1], self.y_centroid-self.semi_lengthY) and near(x[2],
        self.z_centroid-self.semi_lengthZ)):

            y[]