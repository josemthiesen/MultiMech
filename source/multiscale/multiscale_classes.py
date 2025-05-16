# Routine to store classes for different multiscale boundary conditions

from dolfin import *

from abc import ABC, abstractmethod

import numpy as np

import source.tool_box.functional_tools as functional_tools

import source.tool_box.programming_tools as programming_tools

import source.multiscale.multiscale_expressions as multiscale_expressions

########################################################################
#                Boundary conditions's classes templates               #
########################################################################

# Defines a template for the boundary conditions' classes

class BCsClassTemplate(ABC):

    @abstractmethod

    # The following methods have the pass argument only because they 
    # will be defined in the child classes. But they need to be created
    # in the child classes

    def update(self, bilinear_form, linear_form, trial_functionsDict,
    test_functionsDict, boundary_conditions, monolithic_functionSpace):

        pass

# Defines a class for the minimally-constrained boundary condition using 
# first order homogenization

class MinimallyConstrainedFirstOrderBC(BCsClassTemplate):

    def __init__(self,  constrained_fieldName, fields_names, 
    elements_dictionary, mesh_dataClass, macro_quantitiesFilesDict, 
    macro_quantitiesClasses, volume_inverse, fluctuation_field):
        
        # Gets the name of the field to be minimally constrained and its
        # gradient. But removes non ASCII characters and blank spaces so 
        # this names can be used as variables' names

        self.constrained_fieldName = functional_tools.convert_stringToASCII(
        constrained_fieldName)

        self.constrained_gradientFieldName = (constrained_fieldName+"_"+
        "gradient")
        
        # Evaluates the volume of the microscale

        if volume_inverse is None:

            self.volume_inverse = (1.0/assemble(1.0*mesh_dataClass.dx))

        else:

            self.volume_inverse = volume_inverse
        
        # Stores the macro quantities that will be used as boundary con-
        # ditions

        self.macro_quantitiesClasses = macro_quantitiesClasses

        # Swaps the keys of the files dictionary to keep the names of 
        # the field and its gradient as variables in the macro quanti-
        # ties class. Then, updates the list of macro quantities' clas-
        # ses

        self.macro_quantitiesClasses.append(
        functional_tools.MacroQuantitiesInTime(
        programming_tools.change_dictionaryKeys(
        macro_quantitiesFilesDict, [["macro field file", 
        self.constrained_fieldName], ["macro field gradient file",
        self.constrained_gradientFieldName]])))

        # Gets the number of dimensions of the primal field

        n_dimsPrimalField = len(elements_dictionary[
        self.constrained_fieldName].value_shape())

        # Constructs the fields of the Lagrange multipliers to cons-
        # train the field and its gradient

        if n_dimsPrimalField==0:

            self.lagrange_fieldElement = FiniteElement("Real", 
            mesh_dataClass.mesh.ufl_cell(), 0)

            self.lagrange_gradFieldElement = VectorElement("Real", 
            mesh_dataClass.mesh.ufl_cell(), 0)

        elif n_dimsPrimalField==1:

            self.lagrange_fieldElement = VectorElement("Real", 
            mesh_dataClass.mesh.ufl_cell(), 0)

            self.lagrange_gradFieldElement = TensorElement("Real", 
            mesh_dataClass.mesh.ufl_cell(), 0)

        else:

            # Gets the dimension of the space

            space_dimension = mesh_dataClass.mesh.topology().dim()

            # Constructs the tensor element

            dimensions_list = [space_dimension for tensor_index in range(
            n_dimsPrimalField)]

            self.lagrange_fieldElement = TensorElement("Real", 
            mesh_dataClass.mesh.ufl_cell(), 0, tuple(dimensions_list))

            # Adds oen tensorial index due to the gradient

            self.lagrange_gradFieldElement = TensorElement("Real", 
            mesh_dataClass.mesh.ufl_cell(), 0, tuple(
            dimensions_list.append(space_dimension)))
        
        # Gives names for the Lagrange multiplier fiels and adds them to
        # the list of fields' names. These names will be used later to 
        # pick up the trial and test functions relative to these Lagran-
        # ge fields

        self.lagrange_fieldName = (self.constrained_fieldName+"_lagran"+
        "ge_multiplier")

        self.lagrange_gradFieldName = (self.constrained_gradientFieldName
        +"_lagrange_multiplier")

        self.fields_names = fields_names

        self.fields_names.extend([self.lagrange_fieldName, 
        self.lagrange_gradFieldName])
        
        # Adds the elements of the Lagrange multipliers

        self.elements_dictionary = elements_dictionary

        self.elements_dictionary[self.lagrange_fieldName] = (
        self.lagrange_fieldElement)

        self.elements_dictionary[self.lagrange_gradFieldName] = (
        self.lagrange_gradFieldElement)

        # Verifies if fluctuation field flag is true, which means that 
        # the field in the microscale boundary value problem is, in 
        # fact, the fluctuation

        if fluctuation_field:

            self.field_correction = (getattr(
            self.macro_quantitiesClasses[-1], self.constrained_fieldName
            )+dot(getattr(self.macro_quantitiesClasses[-1], 
            self.constrained_gradientFieldName), mesh_dataClass.x))

        else:

            if n_dimsPrimalField==0:

                self.field_correction = Constant(0.0)

            elif n_dimsPrimalField==1:

                self.field_correction = Constant([0.0, 0.0, 0.0])

            else:

                # Gets the dimension of the space

                space_dimension = mesh_dataClass.mesh.topology().dim()

                # Constructs the tensor element

                dimensions_list = [space_dimension for tensor_index in range(
                n_dimsPrimalField)]

                self.field_correction = Constant(np.zeros(tuple(
                dimensions_list)))
        
    # Defines a function to in fact construct the boundary condition

    def update(self, bilinear_form, linear_form, boundary_conditions,
    trial_functionsDict, test_functionsDict, boundary_conditionIndex,
    mesh_dataClass, monolithic_functionSpace):

        # Gets the trial and test functions

        primal_fieldTrial = trial_functionsDict[
        self.constrained_fieldName]

        primal_fieldTest = test_functionsDict[self.constrained_fieldName]

        lagrange_fieldTrial = trial_functionsDict[self.lagrange_fieldName]

        lagrange_fieldTest = test_functionsDict[self.lagrange_fieldName]

        lagrange_gradFieldTrial = trial_functionsDict[
        self.lagrange_gradFieldName]

        lagrange_gradFieldTest = test_functionsDict[
        self.lagrange_gradFieldName]

        # Adds the variational bits to the bilinear form

        # Adds the work done by the constraint on the field

        bilinear_form += (self.volume_inverse*((inner(lagrange_fieldTest,
        getattr(self.macro_quantitiesClasses[boundary_conditionIndex], 
        self.constrained_fieldName))*mesh_dataClass.dx)-(inner(
        lagrange_fieldTest, primal_fieldTrial)*mesh_dataClass.dx)))

        # Adds the work done by the constraint on the gradient of the 
        # field

        bilinear_form += (self.volume_inverse*((inner(
        lagrange_gradFieldTest, getattr(self.macro_quantitiesClasses[
        boundary_conditionIndex], self.constrained_gradientFieldName))*
        mesh_dataClass.dx)-(inner(lagrange_gradFieldTest, grad(
        primal_fieldTrial))*mesh_dataClass.dx)))

        # Adds the work done by the lagrange multipliers

        linear_form += (self.volume_inverse*((inner(lagrange_fieldTrial, 
        primal_fieldTest)*mesh_dataClass.dx)+(inner(grad(primal_fieldTest
        ), lagrange_gradFieldTrial)*mesh_dataClass.dx)))

        # Returns the boundary conditions and variational forms

        return bilinear_form, linear_form, boundary_conditions
    
# Defines a class for the linear boundary condition considering first 
# order homogenization

class LinearFirstOrderBC(BCsClassTemplate):

    def __init__(self,  constrained_fieldName, fields_names, 
    elements_dictionary, mesh_dataClass, macro_quantitiesFilesDict, 
    macro_quantitiesClasses, volume_inverse, fluctuation_field):
        
        # Gets the name of the field to be minimally constrained and its
        # gradient. But removes non ASCII characters and blank spaces so 
        # this names can be used as variables' names

        self.constrained_fieldName = functional_tools.convert_stringToASCII(
        constrained_fieldName)

        self.constrained_gradientFieldName = (constrained_fieldName+"_"+
        "gradient")
        
        # Evaluates the volume of the microscale

        if volume_inverse is None:

            self.volume_inverse = (1.0/assemble(1.0*mesh_dataClass.dx))

        else:

            self.volume_inverse = volume_inverse
        
        # Stores the macro quantities that will be used as boundary con-
        # ditions

        self.macro_quantitiesClasses = macro_quantitiesClasses

        # Swaps the keys of the files dictionary to keep the names of 
        # the field and its gradient as variables in the macro quanti-
        # ties class. Then, updates the list of macro quantities' clas-
        # ses

        self.macro_quantitiesClasses.append(
        functional_tools.MacroQuantitiesInTime(
        programming_tools.change_dictionaryKeys(
        macro_quantitiesFilesDict, [["macro field file", 
        self.constrained_fieldName], ["macro field gradient file",
        self.constrained_gradientFieldName]])))

        # Finds the centroid of the mesh

        position_vector = mesh_dataClass.x

        x_centroid = (self.volume_inverse*assemble(position_vector[
        0]*mesh_dataClass.dx))

        y_centroid = (self.volume_inverse*assemble(position_vector[
        1]*mesh_dataClass.dx))

        z_centroid = (self.volume_inverse*assemble(position_vector[
        2]*mesh_dataClass.dx))

        # Saves the number of fields

        self.n_fields = len(fields_names)

        # Finds the index of this field to be constrained

        self.field_index = 0

        for i in range(self.n_fields):

            if fields_names[i]==self.constrained_fieldName:

                self.field_index = i 

                break

        # Saves the atributes that will be retrieved again for latter
        # classes
        
        self.elements_dictionary = elements_dictionary

        self.fields_names = fields_names

        # Gets the number of dimensions of the primal field

        self.n_dimsPrimalField = len(elements_dictionary[
        self.constrained_fieldName].value_shape())

        # Verifies if fluctuation field flag is true, which means that 
        # the field in the microscale boundary value problem is, in 
        # fact, the fluctuation

        if fluctuation_field:

            self.field_correction = (getattr(
            self.macro_quantitiesClasses[-1], self.constrained_fieldName
            )+dot(getattr(self.macro_quantitiesClasses[-1], 
            self.constrained_gradientFieldName), mesh_dataClass.x))

            # Constructs the expressions for the field on the boundary 
            # considering zero fluctuations on the boundary

            if self.n_dimsPrimalField==0:

                self.field_expression = "0.0"

            elif self.n_dimsPrimalField==1:

                self.field_expression = ("0.0", "0.0", "0.0")

            else:

                raise ValueError("There is not possible yet to use the"+
                " linear multiscale boundary condition for fields with"+
                " dimensionality larger than a 1 (a vector). The given"+
                " dimensionality is "+str(self.n_dimsPrimalField))

        else:

            # Constructs the expressions for the field on the boundary 
            # considering zero fluctuations on the boundary

            if self.n_dimsPrimalField==0:

                self.field_expression = multiscale_expressions.LinearFieldExpression(
                getattr(self.macro_quantitiesClasses[-1], 
                self.constrained_fieldName), getattr(
                self.macro_quantitiesClasses[-1], 
                self.constrained_gradientFieldName), x_centroid, 
                y_centroid, z_centroid)

                # As the field is not the fluctuation, the correction to
                # the full field is null

                self.field_correction = Constant(0.0)

            elif self.n_dimsPrimalField==1:

                self.field_expression = multiscale_expressions.LinearVectorFieldExpression(
                getattr(self.macro_quantitiesClasses[-1], 
                self.constrained_fieldName), getattr(
                self.macro_quantitiesClasses[-1], 
                self.constrained_gradientFieldName), x_centroid, 
                y_centroid, z_centroid)

                # As the field is not the fluctuation, the correction to
                # the full field is null

                self.field_correction = Constant([0.0, 0.0, 0.0])

            else:

                raise ValueError("There is not possible yet to use the"+
                " linear multiscale boundary condition for fields with"+
                " dimensionality larger than a 1 (a vector). The given"+
                " dimensionality is "+str(self.n_dimsPrimalField))
        
    # Defines a function to in fact construct the boundary condition

    def update(self, bilinear_form, linear_form, boundary_conditions,
    trial_functionsDict, test_functionsDict, boundary_conditionIndex, 
    mesh_dataClass, monolithic_functionSpace):
        
        # Adds the zero fluctuation boundary condtion. Checks first if
        # the formulation has multiple fields

        if self.n_fields>1:

            boundary_conditions.append(DirichletBC(
            monolithic_functionSpace.sub(self.field_index), 
            self.field_expression, "on_boundary"))

        else:

            boundary_conditions.append(DirichletBC(
            monolithic_functionSpace, self.field_expression, "on_bound"+
            "ary"))

        # Returns the boundary conditions and variational forms

        return bilinear_form, linear_form, boundary_conditions