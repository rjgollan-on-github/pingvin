#!/usr/bin/python3
"""
Calls the busemann.py module to generate a contour file.

Usage:

  > ./generate-cross-sections.py n_pts_per_section n_sections

.. Author: RJG
.. Date: 2024-08-04
"""

from gdtk.busemann import BusemannDiffuser
from math import asin, pi, sin, cos
import numpy as np
import sys
import os

xsect_dir = "xsect"
if not os.path.exists(xsect_dir):
    os.mkdir(xsect_dir)

if len(sys.argv) != 3:
    print("Usage:")
    print("> ./generate-cross-sections.py n_pts_per_section n_sections")
    sys.exit(1)

n_pts = int(sys.argv[1])
n_sections = int(sys.argv[2])

M2 = 3.0
k = 1.4
theta_23 = asin(k/M2)
bd = BusemannDiffuser(M2, theta_23)

theta_1 = 2.9673531127718715
theta_2 = 0.3026840705911106
dtheta_max = (theta_1 - theta_2)/1000.0
r = 1.0
bd.generate_contour(r, dtheta_max)
print(bd.properties())

bd.write_wall_properties("b-props.dat")

thetas = np.linspace(bd._thetas[-1], bd._thetas[0], n_sections, endpoint=True)
for i, theta in enumerate(thetas):
    u, v, r = bd.ode_soln(theta)
    x, y = bd.xy_from_rtheta(r, theta)
    drdx = bd.dydx(theta, u, v)
    # Build a cross-section
    fname = f"{xsect_dir}/xsect-{i:03d}"
    with open(fname, "w") as f:
        f.write(f"{n_pts}\n")
        dphi = (pi/2)/(n_pts-1)
        for j in range(n_pts):
            phi = j*dphi
            f.write(f"{x:.6e} {y*sin(phi):.6e} {y*cos(phi):.6e} {drdx*sin(phi):.6e} {drdx*cos(phi):.6e}\n")
