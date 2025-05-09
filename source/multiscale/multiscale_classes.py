# Routine to store classes for different multiscale boundary conditions

from dolfin import *

import source.tool_box.functional_tools as functional_tools

########################################################################
#                Boundary conditions's classes templates               #
########################################################################

# Defines a template for the boundary conditions' classes

class BCsClassTemplate:

    def __init__(self, initialization_function, update_function,
    additional_information, code_providedInfo):

        # Each class must have four variables: the initialization func-
        # tion, i.e. the function that initializes all concerning varia-
        # bles to the post-process method; the update function, which 
        # runs the process proper; a list of required additional infor-
        # mation names, which will be used to retrieve information given
        # in dictionaries by the user, thus, the names are their keys; 
        # a list of variables given by the code, that are retrieven from
        # the context class

        self.initialization_function = initialization_function

        self.update_function = update_function

        self.additional_information = additional_information

        self.code_providedInfo = code_providedInfo

# Defines a class for the minimally-constrained boundary condition

@programming_tools.optional_argumentsInitializer({'macro_quantitiesCla'+
'sses': lambda: []})

class MinimallyConstrainedFirstOrderBC:

    def __init__(self, primal_functionSpace, mesh_dataClass, 
    macro_quantitiesFilesDict, macro_quantitiesClasses=None):

        self.macro_quantitiesClasses = macro_quantitiesClasses

        # Updates the list of macro quantities' classes

        self.macro_quantitiesClasses.append(
        functional_tools.MacroQuantitiesInTime(macro_quantitiesFilesDict
        ))

        # Gets the number of dimensions of the primal field

        n_dimsPrimalField = len(primal_functionSpace.ufl_element(
        ).value_shape())

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

            raise ValueError("Minimally constrained BC is not yet impl"+
            "emented for fields that are not scalar or vector fields")