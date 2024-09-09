"""
Calls gmsh python module to create a 2D unstructured quad mesh.

Usage:

    > ./generate-grid-2d.py n_wall_pts bfn

.. Author: RJG
.. Date: 2024-07-07
"""

import gmsh
from math import pi, ceil
import sys

if len(sys.argv) != 3:
    print("Usage:")
    print("> ./generate-grid-2d.py n_wall_pts bfn")
    sys.exit(1)

n_wall_pts = int(sys.argv[1])
bfn = sys.argv[2]

gmsh.initialize()
gmsh.model.add(bfn)

a = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
b = gmsh.model.geo.addPoint(0.0, 1.0, 0.0)
c = gmsh.model.geo.addPoint(1.0, 0.0, 0.0)

ba = gmsh.model.geo.addLine(b, a)
ac = gmsh.model.geo.addLine(a, c)
cb = gmsh.model.geo.addCircleArc(c, a, b)

wall = gmsh.model.addPhysicalGroup(1, [cb])
gmsh.model.setPhysicalName(1, wall, "wall")
symm = gmsh.model.addPhysicalGroup(1, [ba, ac])
gmsh.model.setPhysicalName(1, symm, "symm")

gmsh.option.setNumber("Mesh.SaveAll", 1)

curve = gmsh.model.geo.addCurveLoop([ac, cb, ba])

surf = gmsh.model.geo.addPlaneSurface([curve])

gmsh.model.geo.synchronize()

wall_points = n_wall_pts
symm_points = ceil(wall_points/(pi/2))
print(f"{wall_points=}")
print(f"{symm_points=}")
gmsh.model.mesh.setTransfiniteCurve(ba, symm_points)
gmsh.model.mesh.setTransfiniteCurve(ac, symm_points)
gmsh.model.mesh.setTransfiniteCurve(cb, wall_points)

gmsh.option.setNumber("Mesh.Algorithm", 11)
#gmsh.option.setNumber("Mesh.QuadqsSizemapMethod", 0)
gmsh.option.setNumber("Mesh.QuadqsTopologyOptimizationMethods", 100)
#gmsh.option.setNumber("Mesh.QuadqsRemeshingBoldness", 1.0)
gmsh.option.setNumber("Mesh.QuadqsScalingOnTriangulation", 0.75)
gmsh.model.mesh.generate(2)
gmsh.write(f"{bfn}.su2")
gmsh.write(f"{bfn}.vtk")
