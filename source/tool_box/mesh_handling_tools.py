# Routine to store some methods to work with meshes

from dolfin import *

import copy

########################################################################
#                              Mesh files                              #
########################################################################

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

    data_meshCollection = MeshValueCollection("size_t", mesh, 
    mesh.topology().dim())

    # Reads the mesh with domain physical groups

    with XDMFFile(file_name+"_domain.xdmf") as infile:

        infile.read(mesh)

        infile.read(data_meshCollection, "domain")

    # Converts the mesh value collection to mesh function, for mesh va-
    # lue collections are low level and cannot be used for FEM integra-
    # tion and other higher level operations inside FEniCS

    domain_meshFunction = MeshFunction("size_t", mesh, 
    data_meshCollection)

    # Reinitializes the mesh value colection to the boundary data

    data_meshCollection = MeshValueCollection("size_t", mesh,
    mesh.topology().dim()-1)

    # Reads the mesh with surface physical groups

    with XDMFFile(file_name+"_boundary.xdmf") as infile:
    
        infile.read(data_meshCollection, "boundary")

    # Converts the mesh value collection to mesh function

    boundary_meshFunction = MeshFunction("size_t", mesh, 
    data_meshCollection)

    # Sets the integration differentials

    dx = Measure("dx", domain=mesh, subdomain_data=domain_meshFunction)

    ds = Measure("ds", domain=mesh, subdomain_data=boundary_meshFunction)

    # Sets the normal vector to the mesh's boundary

    n  = FacetNormal(mesh)

    # Returns these objects

    return mesh, dx, ds, n, domain_meshFunction, boundary_meshFunction

########################################################################
#                              Submeshing                              #
########################################################################

# Defines a function to generate a submesh using MeshView and creates

def create_submesh(mesh, original_cellMarkers, volume_physGroupsTags):

    # If the submesh is meant to be constructed from different volume 
    # physical groups of the mesh, the object of cell_markers has to be
    # changed, because the MeshView create function works only with a 
    # single marker to define a subdomain and, from there, a submesh. 
    # Therefore, the physical groups that contain the RVE are collapsed 
    # into a single marker; conventionally the first selected marker

    if isinstance(volume_physGroupsTags, list):

        # Copies the original cell markers object, so that it won't be 
        # affected elsewhere in the code

        cell_markers = copy.deepcopy(original_cellMarkers)

        # Iterates through the elements of the parent mesh

        if len(volume_physGroupsTags)>1:

            for element in cells(mesh):

                # Tests if the cell marker is in the list of required 
                # cell markers

                if cell_markers[element] in volume_physGroupsTags:

                    # Changes the cell marker to the first required one

                    cell_markers[element] = volume_physGroupsTags[0]

        # Collapses the volume physical groups tags list to the first 
        # component

        volume_physGroupsTags = volume_physGroupsTags[0]

    # Creates a submesh for the RVE

    RVE_submesh = MeshView.create(cell_markers,volume_physGroupsTags)

    # Returns the submesh

    return RVE_submesh