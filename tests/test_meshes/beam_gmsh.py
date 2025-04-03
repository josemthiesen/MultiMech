import gmsh

import numpy as np

def generate_beam(length_x, length_y, length_z, transfinite_x, 
transfinite_y, transfinite_z, directory_path, file_name):

    gmsh.initialize()

    # Adds the points

    p11 = gmsh.model.geo.addPoint(length_x, 0.0, 0.0)

    p12 = gmsh.model.geo.addPoint(length_x, length_y, 0.0)

    p13 = gmsh.model.geo.addPoint(0.0, length_y, 0.0)

    p14 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)

    p21 = gmsh.model.geo.addPoint(length_x, 0.0, length_z)

    p22 = gmsh.model.geo.addPoint(length_x, length_y, length_z)

    p23 = gmsh.model.geo.addPoint(0.0, length_y, length_z)

    p24 = gmsh.model.geo.addPoint(0.0, 0.0, length_z)

    gmsh.model.geo.synchronize()

    # Adds the lines

    l11 = gmsh.model.geo.addLine(p11, p12)

    l12 = gmsh.model.geo.addLine(p12, p13)

    l13 = gmsh.model.geo.addLine(p13, p14)

    l14 = gmsh.model.geo.addLine(p14, p11)

    l21 = gmsh.model.geo.addLine(p21, p22)

    l22 = gmsh.model.geo.addLine(p22, p23)

    l23 = gmsh.model.geo.addLine(p23, p24)

    l24 = gmsh.model.geo.addLine(p24, p21)

    l121 = gmsh.model.geo.addLine(p11, p21)

    l122 = gmsh.model.geo.addLine(p12, p22)

    l123 = gmsh.model.geo.addLine(p13, p23)

    l124 = gmsh.model.geo.addLine(p14, p24)

    gmsh.model.geo.synchronize()

    # Adds the surfaces

    surface_1 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop([l11, l12, l13, l14])])

    surface_2 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop([l11, l122, -l21, -l121])])

    surface_3 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop([l12, l123, -l22, -l122])])

    surface_4 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop([l13, l124, -l23, -l123])])

    surface_5 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop([l14, l121, -l24, -l124])])

    surface_6 = gmsh.model.geo.addPlaneSurface([
    gmsh.model.geo.addCurveLoop([l21, l22, l23, l24])])

    gmsh.model.geo.synchronize()

    # Adds the volume

    volume = gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop([
    surface_1, surface_2, surface_3, surface_4, surface_5, surface_6])])

    gmsh.model.geo.synchronize()

    # Makes the lines transfinite

    gmsh.model.geo.mesh.setTransfiniteCurve(l11, transfinite_y)

    gmsh.model.geo.mesh.setTransfiniteCurve(l21, transfinite_y)

    gmsh.model.geo.mesh.setTransfiniteCurve(l13, transfinite_y)

    gmsh.model.geo.mesh.setTransfiniteCurve(l23, transfinite_y)

    gmsh.model.geo.mesh.setTransfiniteCurve(l12, transfinite_x)

    gmsh.model.geo.mesh.setTransfiniteCurve(l22, transfinite_x)

    gmsh.model.geo.mesh.setTransfiniteCurve(l14, transfinite_x)

    gmsh.model.geo.mesh.setTransfiniteCurve(l24, transfinite_x)

    gmsh.model.geo.mesh.setTransfiniteCurve(l121, transfinite_z)

    gmsh.model.geo.mesh.setTransfiniteCurve(l122, transfinite_z)

    gmsh.model.geo.mesh.setTransfiniteCurve(l123, transfinite_z)

    gmsh.model.geo.mesh.setTransfiniteCurve(l124, transfinite_z)

    gmsh.model.geo.synchronize()

    # Makes the surfaces transfinite

    gmsh.model.geo.mesh.setTransfiniteSurface(surface_1)

    gmsh.model.geo.mesh.setTransfiniteSurface(surface_2)

    gmsh.model.geo.mesh.setTransfiniteSurface(surface_3)

    gmsh.model.geo.mesh.setTransfiniteSurface(surface_4)

    gmsh.model.geo.mesh.setTransfiniteSurface(surface_5)

    gmsh.model.geo.mesh.setTransfiniteSurface(surface_6)

    gmsh.model.geo.synchronize()

    # Makes the volume transfinite

    gmsh.model.geo.mesh.setTransfiniteVolume(volume)

    gmsh.model.geo.synchronize()
    
    gmsh.model.geo.removeAllDuplicates()

    gmsh.model.geo.synchronize()

    # Adds the physical groups
    
    gmsh.model.addPhysicalGroup(3, volume, 1, name='Generic volume')
    
    gmsh.model.addPhysicalGroup(2, surface_1, 2, name='back')
    
    gmsh.model.addPhysicalGroup(2, surface_2, 3, name='right')
    
    gmsh.model.addPhysicalGroup(2, surface_3, 4, name='lower')
    
    gmsh.model.addPhysicalGroup(2, surface_4, 5, name='left')
    
    gmsh.model.addPhysicalGroup(2, surface_5, 6, name='upper')
    
    gmsh.model.addPhysicalGroup(2, surface_6, 7, name='front')

    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(3)

    gmsh.write(directory_path+"//"+file_name+".msh")
    
    gmsh.fltk.run()

    gmsh.finalize()

# Defines a function to generate the beam given the micropolar model's
# parameters

def generate_micropolarBeam(mu, ratio_Lb, beta, gamma, directory_path,
file_name):

    # Evaluate the characteristic length and the corresponding beam's 
    # square cross section dimension

    L = np.sqrt((beta+gamma)/(2*mu))

    b = L/ratio_Lb

    # Defines the mesh's number of divions

    transfinite_x = 7

    transfinite_y = 7
    
    transfinite_z = 21

    # Generates the mesh

    generate_beam(b, b, 10*b, transfinite_x, transfinite_y, 
    transfinite_z, directory_path, file_name)