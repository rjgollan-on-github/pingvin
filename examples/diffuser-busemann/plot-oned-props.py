# Author: RJG
# Date: 2024-09-23
#
# To use:
#   > python plot-oned-props.py

import numpy as np
import matplotlib.pyplot as plt

FNAME = "pingvin-oned-stream-thrust-avg.data"

# Only interested in a subset of data, so pick that out by column
x, area, p, Mach = np.loadtxt(FNAME, skiprows=1, usecols=(0, 1, 2, 5), unpack=True)

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.75)
twin1 = ax.twinx()
twin2 = ax.twinx()
twin2.spines['right'].set_position(("axes", 1.125))

p0, = ax.plot(x, p/1000.0, label='pressure, kPa', color='red')
p1, = twin1.plot(x, Mach, label='Mach number', color='blue')
p2, = twin2.plot(x, area, label='area, m^2', color='black')

twin1.yaxis.label.set_color(p1.get_color())
twin2.yaxis.label.set_color(p2.get_color())

ax.legend(handles=[p0, p1, p2], loc='center left', frameon=True, framealpha=1, edgecolor="black")
ax.set_ylabel("pressure, kPa")
ax.set_xlabel("x, m")
twin1.set_ylabel("Mach number")
twin2.set_ylabel("area, m^2")
twin1.set_ylim(0, 6)
twin2.set_ylim(0, 1.0)
tkw = dict(size=4, width=1.5)
twin1.tick_params(axis='y', colors=p1.get_color(), **tkw)
twin2.tick_params(axis='y', colors=p2.get_color(), **tkw)


plt.savefig("busemann-oned-props.png", dpi=600)

