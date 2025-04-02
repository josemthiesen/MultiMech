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

    # Recovers the physical groups

    domain_physicalGroupsNameToTag = dict()

    boundary_physicalGroupsNameToTag = dict()

    # Gets the physical group through the field data

    field_data = mesh_reading.field_data

    # Iterates through its keys

    for physical_group, tag_dimensionality in field_data.items():

        # If the dimensionality is 2, gets the boundary

        if tag_dimensionality[1]==2:

            boundary_physicalGroupsNameToTag[physical_group] = tag_dimensionality[
            0]

        # If the dimensionality is 3, gets the domain

        elif tag_dimensionality[1]==3:

            domain_physicalGroupsNameToTag[physical_group] = tag_dimensionality[
            0]

    print("###########################################################"+
    "#############\n#                        Mesh - physical groups   "+
    "                     #\n#########################################"+
    "###############################\n")

    print("Finds the following domain physical groups with their respe"+
    "ctive tags:")

    for physical_group, tag in domain_physicalGroupsNameToTag.items():

        print(physical_group, ": ", tag)

    print("\nFinds the following boundary physical groups with their r"+
    "espective tags:")

    for physical_group, tag in boundary_physicalGroupsNameToTag.items():

        print(physical_group, ": ", tag)

    print("")

    # Gets the cells which consist of the desired element

    for i in range(len(desired_elements)):

        print("Saves the mesh of", data_sets[i], "dataset\n")

        cells_dictionary[desired_elements[i]] = (
        mesh_reading.get_cells_type(desired_elements[i]))

        print("There are "+str(len(cells_dictionary[desired_elements[i]]
        ))+" "+str(desired_elements[i])+" elements in the mesh.\n")

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

    print("###########################################################"+
    "#############\n#                    Mesh reading has been finaliz"+
    "ed                   #\n#########################################"+
    "###############################\n")

    # Then, calls the xdmf reader and returns its output

    return read_xdmfMesh(file_name, domain_physicalGroupsNameToTag=
    domain_physicalGroupsNameToTag, boundary_physicalGroupsNameToTag=
    boundary_physicalGroupsNameToTag)

# Defines a function to read a mesh from a xdmf file

def read_xdmfMesh(file_name, domain_physicalGroupsNameToTag=dict(), 
boundary_physicalGroupsNameToTag=dict()):
    
    if len(domain_physicalGroupsNameToTag.keys())==0 or (len(
    boundary_physicalGroupsNameToTag.keys())==0):
        
        print("WARNING: the dictionaries of physical groups' names to "+
        "tags are empty, thus, it won't be possible to use the names i"+
        "n the variational forms. Use a .msh mesh file directly instea"+
        "d.")

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
    boundary_meshCollection, boundary_meshFunction, 
    domain_physicalGroupsNameToTag, boundary_physicalGroupsNameToTag)

########################################################################
#                              Submeshing                              #
########################################################################

# Defines a function to generate a submesh using MeshView and creates

def create_submesh(domain_meshCollection, volume_physGroupsTags, 
parent_functionSpace):

    # Gets the mesh, the polynomial degree, and the shape function

    mesh = parent_functionSpace.mesh()

    polynomial_degree = parent_functionSpace.ufl_element().degree()

    shape_function = parent_functionSpace.ufl_element().family()

    # Gets the dimensionality of the parent function space

    n_dimsParentFunctionSpace = len(parent_functionSpace.ufl_element(
    ).value_shape())

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

    # If the element is mixed, the .family() will return 'Mixed'

    if shape_function=='Mixed':

        submesh_functionSpace = FunctionSpace(sub_mesh,
        parent_functionSpace.ufl_element())

    else:

        # If the field is scalar, the dimensionality is 0

        if n_dimsParentFunctionSpace==0:

            submesh_functionSpace = FunctionSpace(sub_mesh, 
            shape_function, polynomial_degree)

        # If the field is vector function, the dimensionality is 1

        elif n_dimsParentFunctionSpace==1:

            submesh_functionSpace = VectorFunctionSpace(sub_mesh,
            shape_function, polynomial_degree)

        # If the field is a second order tensor function

        elif n_dimsParentFunctionSpace==2:

            submesh_functionSpace = TensorFunctionSpace(sub_mesh, 
            shape_function, polynomial_degree)

        # Handle not implemented function spaces

        else:

            raise NameError("create_submesh in mesh_handling_tools.py "+
            "does not support a function space with dimensionality lar"+
            "ger than "+str(n_dimsParentFunctionSpace)+", for it's not"+
            " been implemented yet or it's not possible to implement.")
        
    # Initializes the function to the solution at the submesh

    submesh_function = Function(submesh_functionSpace)

    # Initializes the DOF mappings for the RVE and for the original mesh

    sub_meshMapping = []

    parent_meshMapping = []

    # Verifies whether there is only on field

    if shape_function!='Mixed':

        sub_meshMapping.append(submesh_functionSpace.dofmap())

        parent_meshMapping.append(parent_functionSpace.dofmap())

    # If there are multiple fields

    else:

        # Iterates through the number of fields

        for i in range(parent_functionSpace.ufl_element().num_sub_elements()):

            # Adds the submesh mapping and the parent mesh mapping

            sub_meshMapping.append(submesh_functionSpace.sub(i).dofmap())

            parent_meshMapping.append(parent_functionSpace.sub(i).dofmap(
            ))
            
    # Sets the integration differential in the submesh

    dx_submesh = Measure("dx", domain=sub_mesh)

    # Returns the submesh, the updated cell markers, and the DOF mappings

    return (sub_mesh, submesh_cellMarkers, submesh_functionSpace, 
    sub_meshMapping, parent_meshMapping, submesh_function, 
    sub_toParentCellMap, dx_submesh)

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
