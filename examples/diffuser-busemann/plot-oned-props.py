# Author: RJG
# Date: 2024-09-23
#
# To use:
#   > python plot-oned-props.py

import numpy as np
import matplotlib.pyplot as plt

FNAME = "pingvin-oned-cmes-avg.data"

# Only interested in a subset of data, so pick that out by column
x, area, p, Mach, pt = np.loadtxt(FNAME, skiprows=1, usecols=(0, 1, 2, 5, 7), unpack=True)

fig, ax = plt.subplots()
fig.subplots_adjust(right=0.75)
twin1 = ax.twinx()
twin2 = ax.twinx()
twin2.spines['right'].set_position(("axes", 1.125))

p0, = ax.plot(x, p/1000.0, label='pressure, kPa', color='red')
p1, = twin1.plot(x, Mach, label='Mach number', color='blue')
p2, = ax.plot(x, area/area[-1], label='contraction ratio', color='black')
p3, = twin2.plot(x, pt/pt[0], label='total pressure recovery', color='green')

twin1.yaxis.label.set_color(p1.get_color())
twin2.yaxis.label.set_color(p3.get_color())

ax.legend(handles=[p0, p1, p2, p3], loc=(0.05,0.7), frameon=True, framealpha=1, edgecolor="black", prop={'size':6})
ax.set_ylabel("pressure, kPa; contraction ratio")
ax.set_xlabel("x, m")
twin1.set_ylabel("Mach number")
twin2.set_ylabel("total pressure recovery")
twin1.set_ylim(0, 6)
twin2.set_ylim(0, 1.1)
tkw = dict(size=4, width=1.5)
twin1.tick_params(axis='y', colors=p1.get_color(), **tkw)
twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)


plt.savefig("busemann-oned-props.png", dpi=600)

