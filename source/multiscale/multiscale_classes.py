# Routine to store classes for different multiscale boundary conditions

from dolfin import *

from abc import ABC, abstractmethod

import source.tool_box.functional_tools as functional_tools

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
    test_functionsDict, boundary_conditions):

        pass

# Defines a class for the minimally-constrained boundary condition

class MinimallyConstrainedFirstOrderBC(BCsClassTemplate):

    def __init__(self,  constrained_fieldName, 
    constrained_gradientFieldName, fields_names, elements_dictionary, 
    mesh_dataClass, macro_quantitiesFilesDict, macro_quantitiesClasses):
        
        # Gets the name of the field to be minimally constrained and its
        # gradient

        self.constrained_fieldName = constrained_fieldName

        self.constrained_gradientFieldName = (
        constrained_gradientFieldName)
        
        # Evaluates the volume of the microscale

        self.volume_inverse = (1.0/assemble(1.0*mesh_dataClass.dx))
        
        # Stores the macro quantities that will be used as boundary con-
        # ditions

        self.macro_quantitiesClasses = macro_quantitiesClasses

        # Stores the dictionary of macro quantities

        self.macro_quantitiesFilesDict = macro_quantitiesFilesDict

        # Updates the list of macro quantities' classes

        self.macro_quantitiesClasses.append(
        functional_tools.MacroQuantitiesInTime(macro_quantitiesFilesDict
        ))

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

        self.lagrange_fieldName = (self.constrained_fieldName+" lagran"+
        "ge multiplier")

        self.lagrange_gradFieldName = (self.constrained_gradientFieldName
        +" lagrange multiplier")

        self.fields_names = fields_names

        self.fields_names.extend([self.lagrange_fieldName, 
        self.lagrange_gradFieldName])
        
        # Adds the elements of the Lagrange multipliers

        self.elements_dictionary = elements_dictionary

        self.elements_dictionary[self.lagrange_fieldName] = (
        self.lagrange_fieldElement)

        self.elements_dictionary[self.lagrange_gradFieldName] = (
        self.lagrange_gradFieldElement)
        
    # Defines a function to in fact construct the boundary condition

    def update(self, bilinear_form, linear_form, boundary_conditions,
    trial_functionsDict, test_functionsDict, boundary_conditionIndex,
    mesh_dataClass):

        # Gets the trial and test functions

        primal_fieldTrial = trial_functionsDict[
        self.constrained_fieldName]

        primal_fieldTest = test_functionsDict[self.constrained_fieldName]

        lagrange_fieldTrial = trial_functionsDict[self.lagrange_fieldName]

        lagrange_fieldTest = test_functionsDict[self.lagrange_fieldName]

        lagrange_gradFieldTrial = trial_functionsDict[
        self.constrained_fieldName]

        lagrange_gradFieldTest = test_functionsDict[
        self.constrained_fieldName]

        # Adds the variational bits to the bilinear form

        # Adds the work done by the constraint on the field

        bilinear_form += (self.volume_inverse*((inner(lagrange_fieldTest,
        getattr(self.macro_quantitiesClasses[boundary_conditionIndex], 
        functional_tools.convert_stringToASCII(self.constrained_fieldName
        )))*mesh_dataClass.dx)-(inner(lagrange_gradFieldTest, 
        primal_fieldTrial)*mesh_dataClass.dx)))

        # Adds the work done by the constraint on the gradient of the 
        # field

        bilinear_form += (self.volume_inverse*((inner(
        lagrange_gradFieldTest, getattr(
        functional_tools.convert_stringToASCII(
        self.macro_quantitiesClasses[boundary_conditionIndex], 
        self.constrained_gradientFieldName)))*mesh_dataClass.dx)-(inner(
        lagrange_gradFieldTest, grad(primal_fieldTrial))*mesh_dataClass.dx)))

        # Adds the work done by the lagrange multipliers

        linear_form += (self.volume_inverse*((inner(lagrange_fieldTrial, 
        primal_fieldTest)*mesh_dataClass.dx)+(inner(grad(primal_fieldTest
        ), lagrange_gradFieldTrial)*mesh_dataClass.dx)))

        # Returns the boundary conditions and variational forms

        return bilinear_form, linear_form, boundary_conditions