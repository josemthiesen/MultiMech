# Routine to store user expressions for multiscale boundary conditions

from dolfin import *

import numpy as np

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

        value = (self.mean_field+(self.gradient_field[0]*(x[0]-
        self.x_centroid))+(self.gradient_field[1]*(x[1]-
        self.y_centroid))+(self.gradient_field[2]*(x[2]-
        self.z_centroid)))

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

        self.centroid = Constant((x_centroid, y_centroid, z_centroid))

        print(x_centroid, y_centroid, z_centroid)

    # Defines a function to evaluate the expression

    def eval(self, value, x):

        # Transforms the position vector into a numpy array and shifts
        # it by the centroid

        """value[0] = (self.mean_field[0]+(self.gradient_field[0,0]*x[0])+
        (self.gradient_field[0,1]*x[1])+(self.gradient_field[0,2]*x[2]
        ))

        value[1] = (self.mean_field[1]+(self.gradient_field[1,0]*x[0])+
        (self.gradient_field[1,1]*x[1])+(self.gradient_field[1,2]*x[2]
        ))

        value[2] = (self.mean_field[2]+(self.gradient_field[2,0]*x[0])+
        (self.gradient_field[2,1]*x[1])+(self.gradient_field[2,2]*x[2]
        ))"""

        value = self.mean_field+dot(self.gradient_field,as_vector(x)-
        self.centroid)

    def value_shape(self):

        return (3,)

# Defines a class for an expression for a generic linear field with res-
# pect to space coordinates

class LinearFieldExpression(UserExpression):

    field_shape = (3,)

    def __init__(self, mean_field, gradient_field, x_centroid, 
    y_centroid, z_centroid):

        super().__init__()

        # Initializes the values

        self.mean_field = np.asarray(mean_field, dtype=np.float64)

        self.gradient_field = np.asarray(gradient_field, dtype=
        np.float64)

        self.centroid = np.asarray([x_centroid, y_centroid, z_centroid], 
        dtype=np.float64)

        # Automatically determines the shape of the field

        self.field_shape = self.mean_field.shape

    # Defines a function to evaluate the expression

    def eval(self, value, x):

        # Transforms the position vector into a numpy array and shifts
        # it by the centroid

        x_shifted = np.asarray(x)-self.centroid

        # Evaluates the linear function for a scalar field

        if self.field_shape==():

            value = self.mean_field+np.dot(self.gradient_field, 
            x_shifted)

        # If the field is a vector or a tensor field

        else:

            value[:] = self.mean_field+np.tensordot(self.gradient_field,
            x_shifted, axes=1)

    def value_shape(self):

        return self.field_shape