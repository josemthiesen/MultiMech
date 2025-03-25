# Routine to store some methods to work with meshes

from dolfin import *

import meshio

########################################################################
#                              Mesh files                              #
########################################################################

# Defines a function to read a mesh from a msh file

def read_mshMesh(file_name, desired_elements=['tetra', 'triangle'],
data_sets=["domain", "boundary"]):

    # Reads the saved gmsh mesh using meshio

    mesh_reading = meshio.read(file_name+".msh")

    # Initializes the dictionary of cells and the list of cell data

    cells_dictionary = dict()

    cell_dataList = []

    # Gets the cells which consist of the desired element

    for i in range(len(desired_elements)):

        print("Saves the mesh of", data_sets[i], "dataset\n")

        cells_dictionary[desired_elements[i]] = (
        mesh_reading.get_cells_type(desired_elements[i]))

        # Gets the physical cell data for this element

        cell_dataPhysical = mesh_reading.get_cell_data("gmsh:physical", 
        desired_elements[i])

        # Gets the geometric cell data for this element

        cell_dataGeometric = mesh_reading.get_cell_data("gmsh:geometri"+
        "cal", desired_elements[i])

        # Adds the geometric cell data to the list of cell data

        cell_dataList.append(cell_dataGeometric)

        # Creates a mesh for this data set and saves it

        mesh_set = meshio.Mesh(points=mesh_reading.points, cells={
        desired_elements[i]: cells_dictionary[desired_elements[i]]},
        cell_data={data_sets[i]: [cell_dataPhysical]})

        meshio.write(file_name+"_"+data_sets[i]+".xdmf", mesh_set)

    # Creates the mesh with all information, including geometric infor-
    # mation and saves it
    
    whole_mesh = meshio.Mesh(points=mesh_reading.points, cells=
    cells_dictionary, cell_data={"whole_mesh": cell_dataList})

    meshio.write(file_name+".xdmf", whole_mesh)

    # Then, calls the xdmf reader and returns its output

    return read_xdmfMesh(file_name)

# Defines a function to read a mesh from a xdmf file

def read_xdmfMesh(file_name):

    # Verifies whether there is an extension in the file name

    if ("." in file_name) or ('.' in file_name):

        raise NameError("The name of the mesh file, "+file_name+", has"+
        " an extension. This name musn't have a extension, for xdmf fi"+
        "les only can be read with this method.\n")

    # Initializes the mesh object and reads the xdmf file

    mesh = Mesh()

    # Initializes a mesh value collection to store mesh data of the do-
    # main

    domain_meshCollection = MeshValueCollection("size_t", mesh, 
    mesh.topology().dim())

    # Reads the mesh with domain physical groups

    with XDMFFile(file_name+"_domain.xdmf") as infile:

        infile.read(mesh)

        infile.read(domain_meshCollection, "domain")

    # Converts the mesh value collection to mesh function, for mesh va-
    # lue collections are low level and cannot be used for FEM integra-
    # tion and other higher level operations inside FEniCS

    domain_meshFunction = MeshFunction("size_t", mesh, 
    domain_meshCollection)

    # Reinitializes the mesh value colection to the boundary data

    boundary_meshCollection = MeshValueCollection("size_t", mesh,
    mesh.topology().dim()-1)

    # Reads the mesh with surface physical groups

    with XDMFFile(file_name+"_boundary.xdmf") as infile:
    
        infile.read(boundary_meshCollection, "boundary")

    # Converts the mesh value collection to mesh function

    boundary_meshFunction = MeshFunction("size_t", mesh, 
    boundary_meshCollection)

    # Sets the integration differentials

    dx = Measure("dx", domain=mesh, subdomain_data=domain_meshFunction)

    ds = Measure("ds", domain=mesh, subdomain_data=boundary_meshFunction)

    # Sets the normal vector to the mesh's boundary

    n  = FacetNormal(mesh)

    # Returns these objects

    return (mesh, dx, ds, n, domain_meshCollection, domain_meshFunction, 
    boundary_meshCollection, boundary_meshFunction)

########################################################################
#                              Submeshing                              #
########################################################################

# Defines a function to generate a submesh using MeshView and creates

def create_submesh(mesh, domain_meshCollection, volume_physGroupsTags, 
parent_functionSpace, mixed_element=None, polynomial_degree=None, 
function_spaceType=None):

    # Creates a new cell markers object to not affect the mesh function
    # at other parts of the code

    submesh_cellMarkers = cpp.mesh.MeshFunctionSizet(mesh, 
    domain_meshCollection)

    # If the submesh is meant to be constructed from different volume 
    # physical groups of the mesh, the object of cell_markers has to be
    # changed, because the MeshView create function works only with a 
    # single marker to define a subdomain and, from there, a submesh. 
    # Therefore, the physical groups that contain the RVE are collapsed 
    # into a single marker; conventionally the first selected marker

    if isinstance(volume_physGroupsTags, list):

        # Iterates through the elements of the parent mesh

        if len(volume_physGroupsTags)>1:

            for element in cells(mesh):

                # Tests if the cell marker is in the list of required 
                # cell markers

                if submesh_cellMarkers[element] in volume_physGroupsTags:

                    # Changes the cell marker to the first required one

                    submesh_cellMarkers[element] = volume_physGroupsTags[
                    0]

        # Collapses the volume physical groups tags list to the first 
        # component

        volume_physGroupsTags = volume_physGroupsTags[0]

    # Creates a submesh for the RVE

    sub_mesh = MeshView.create(submesh_cellMarkers, volume_physGroupsTags)

    # Creates the mapping of elements from the submesh to the parent 
    # mesh, i.e. given the element index in the submesh, it throws the 
    # index in the parent mesh

    sub_toParentCellMap = sub_mesh.topology().mapping()[mesh.id()
    ].cell_map()

    # Creates the function spaces

    submesh_functionSpace = 0

    # If the mixed element is not None, it is used

    if mixed_element!=None:

        submesh_functionSpace = FunctionSpace(sub_mesh,mixed_element)

    else:

        # Handles the exception of the polynomial degree being not given

        if (polynomial_degree==None or (not isinstance(polynomial_degree, 
        int))):
            
            raise ValueError("In creat_submesh, mesh_handling_tools.py"+
            ", if mixed_element was not prescribed, the function space"+
            " is supposed as 'CG' per default. Thus, a polynomial degr"+
            "ee must be provided, which wasn't: polynomial_degree="+str(
            polynomial_degree))

        # If the field is scalar

        if function_spaceType=="scalar":

            submesh_functionSpace = FunctionSpace(sub_mesh, "CG",
            polynomial_degree)

        # If the field is vector function

        elif function_spaceType=="vector":

            submesh_functionSpace = VectorFunctionSpace(sub_mesh, "CG", 
            polynomial_degree)

        # If the field is a tensor function

        elif function_spaceType=="tensor":

            submesh_functionSpace = TensorFunctionSpace(sub_mesh, "CG", 
            polynomial_degree)

        # Handle not implemented function spaces

        else:

            raise NameError("create_submesh in mesh_handling_tools.py "+
            "does not support "+str(function_spaceType)+" function spa"+
            "ce, for it's not been implemented yet or it's not possibl"+
            "e to implement.")
        
    # Initializes the function to the solution at the submesh

    submesh_function = Function(submesh_functionSpace)

    # Initializes the DOF mappings for the RVE and for the original mesh

    sub_meshMapping = []

    parent_meshMapping = []

    # Verifies whether there is only on field

    if mixed_element==None:

        sub_meshMapping.append(submesh_functionSpace.dofmap())

        parent_meshMapping.append(parent_functionSpace.dofmap())

    # If there are multiple fields

    else:

        # Iterates through the number of fields

        for i in range(mixed_element.num_sub_elements()):

            # Adds the submesh mapping and the parent mesh mapping

            sub_meshMapping.append(submesh_functionSpace.sub(i).dofmap())

            parent_meshMapping.append(parent_functionSpace.sub(i).dofmap(
            ))

    # Returns the submesh, the updated cell markers, and the DOF mappings

    return (sub_mesh, submesh_cellMarkers, submesh_functionSpace, 
    sub_meshMapping, parent_meshMapping, submesh_function, 
    sub_toParentCellMap)

# Defines a function to update the field parameters vector of a submesh 
# given the corresponding vector at the parent mesh

def field_parentToSubmesh(submesh, field_submesh, field_parentMesh, 
sub_toParentCellMap, sub_meshMapping, parent_meshMapping):

    # Iterates through the elements of the submesh

    for element in cells(submesh):

        # Gets the index of the element in the submesh

        submesh_index = element.index()

        # Gets the index of the element in the parent mesh

        parent_index = sub_toParentCellMap[submesh_index]

        # Iterates through the fields of the meshes. If the mesh has on-
        # ly one field, displacement for instance, the list of DOF map-
        # ping has only one component

        for i in range(len(sub_meshMapping)):

            # Translates the values of the solution using the DOFs map-
            # ping

            field_submesh.vector()[sub_meshMapping[i].cell_dofs(
            submesh_index)] = field_parentMesh.vector()[
            parent_meshMapping[i].cell_dofs(parent_index)] 

    # Returns the submesh field

    return field_submesh
