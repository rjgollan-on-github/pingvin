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

r = 1.0
dtheta = 0.001
bd.generate_contour(r, dtheta)

y_inlet = bd._ys[-1]
scale = 1.0/y_inlet
x_shift = -bd._xs[-1]*scale
x_start = bd._xs[-1]*scale + x_shift
x_end = bd._xs[0]*scale + x_shift
xs = np.linspace(x_start, x_end, n_sections, endpoint=True)
for ix, x in enumerate(xs):
    # find x value and weight
    w = 0
    i_save = 0
    for i in range(len(bd._xs)-1, 0, -1):
        x0 = bd._xs[i]*scale + x_shift
        x1 = bd._xs[i-1]*scale + x_shift
        if x0 <= x <= x1:
            w = (x-x0)/(x1-x0)
            i_save = i
    # Build a cross-section
    fname = f"{xsect_dir}/xsect-{ix:03d}"
    with open(fname, "w") as f:
        f.write(f"{n_pts}\n")
        dtheta = (pi/2)/(n_pts-1)
        for j in range(n_pts):
            theta = j*dtheta
            r = (1 - w)*bd._ys[i_save] + w*bd._ys[i_save-1]
            r *= scale
            drdx = (1 - w)*bd._dydxs[i_save] + w*bd._dydxs[i_save-1]
            f.write(f"{x:.6e} {r*sin(theta):.6e} {r*cos(theta):.6e} {drdx*sin(theta):.6e} {drdx*cos(theta):.6e}\n")
