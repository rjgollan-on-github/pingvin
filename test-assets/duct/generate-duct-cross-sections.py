#!/usr/bin/python3
"""
Calls the busemann.py module to generate a contour file.

Usage:

  > ./generate-duct-cross-sections.py npts_per_xsect n_xsect

.. Author: RJG
.. Date: 2024-09-01
"""

R = 1.0
L = 1.0

from math import asin, pi, sin, cos
import numpy as np
import sys
import os

xsect_dir = "xsect"
if not os.path.exists(xsect_dir):
    os.mkdir(xsect_dir)

if len(sys.argv) != 3:
    print("Usage:")
    print("> ./generate-duct-cross-sections.py n_pts_per_section n_sections")
    sys.exit(1)

n_pts = int(sys.argv[1])
n_sections = int(sys.argv[2])

xs = np.linspace(0.0, L, n_sections, endpoint=True)
for ix, x in enumerate(xs):
    # Build a cross-section
    fname = f"{xsect_dir}/xsect-{ix:03d}"
    with open(fname, "w") as f:
        f.write(f"{n_pts}\n")
        dtheta = (pi/2)/(n_pts-1)
        for j in range(n_pts):
            theta = j*dtheta
            f.write(f"{x:.6e} {R*sin(theta):.6e} {R*cos(theta):.6e} {0.0:.6e} {0.0:.6e}\n")
