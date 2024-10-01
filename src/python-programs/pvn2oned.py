# Author: RJG
# Date: 2024-09-22
#
# History:
#   2024-10-01 -- Change from conserved mass-momentum-energy to
#                 conserved mass-energy-entropy
#

import click
import pyvista as pv
from numpy import median, mean
import sys
from math import sqrt, log
from scipy.optimize import root_scalar

VTUFILE = "pvnsim/flowfield.vtu"
GAMMA = 1.4
R = 287.0
T_ref = 298.15 # K
p_ref = 101.325e3 # Pa

Fields = ['xvel', 'yvel', 'zvel', 'Mach',
          'p', 'T', 'rho']

@click.command()
@click.option("-@", "--at-x-loc", "xLoc",
              help="x location for extraction of single set of one-d properties",
              show_default=False)
@click.option("-n", "--n-slices", "nSlices",
              default=-1,
              help="number of slices for extracting one-d properties (distributed equally in x)",
              show_default=False)
@click.option("-o", "--output-file", "outfile",
              default='pingvin-oned',
              help="base filename for output of slice data",
              show_default=True)
def run(xLoc, nSlices, outfile):

    # Determine what actions to take
    printSingleLoc = False
    doSlicing = False
    if xLoc != None:
        xLoc = float(xLoc)
        printSingleLoc = True
    if nSlices >= 1:
        doSlicing = True

    # Load data file (or just exit)
    if doSlicing or printSingleLoc:
        mesh = pv.read(VTUFILE)
    else:
        print("pvn2oned: No action requested.")
        print("pvn2oned: Consider using one or both of:")
        print("   --at-x-loc X.Y")
        print("   --n-slices N")
        print("pvnoned: Exiting.")
        return

    if printSingleLoc:
        slice = mesh.slice(origin=(xLoc, 0.0, 0.0), normal=(1.0, 0.0, 0.0))
        # Compute GAMMA and R from information slice
        set_gamma_and_r(slice)
        area_avg, mass_flow_avg, cmes_avg, area = compute_slice_averages(slice, xLoc)

        
        print(f"# Averaged properties at x = {xLoc:.3f}")
        print("---")
        print(f"x-location:      [{xLoc:.6e}, 'm']")
        print(f"area-slice:      [{area:.6e}, 'm^2']")
        print(f"number-of-cells:  {slice.n_cells}")
        print("")
        print("conserved-mass-energy-entropy-average:")
        print_avg(cmes_avg)
        print("")
        print("mass-flow-weighted-average:")
        print_avg(mass_flow_avg)
        print("")
        print("area-weighted-average:")
        print_avg(area_avg)
        print("...")

    if doSlicing:
        slices = mesh.slice_along_axis(n=nSlices, axis="x", tolerance=0.001)
        vars = ('p', 'T', 'rho', 'Mach', 'xvel', 'p0', 'T0', 'rho0')

        with open(f"{outfile}-area-avg.data", "w") as fa, \
             open(f"{outfile}-mass-flow-avg.data", "w") as fm, \
             open(f"{outfile}-cmes-avg.data", "w") as fs:

            header = "x area "
            for v in vars: header += v + " "
            header += "\n"
            fa.write(header)
            fm.write(header)
            fs.write(header)
             
            for i, slice in enumerate(slices):
                if i == 0:
                    set_gamma_and_r(slice)
                x = slice.get_cell(0).center[0]
                print(f"-- slice-{i:03d} @ x={x:.3f} --")
                print("     computing averaged properties")
                area_avg, mass_flow_avg, cmes_avg, area = compute_slice_averages(slice, xLoc)
                x_area_str = f"{x:.6e} {area:.6e} "
                fa.write(x_area_str)
                fm.write(x_area_str)
                fs.write(x_area_str)
                for v in vars:
                    fa.write(f"{area_avg[v]:.6e} ")
                    fm.write(f"{mass_flow_avg[v]:.6e} ")
                    fs.write(f"{cmes_avg[v]:.6e} ")
                fa.write("\n")
                fm.write("\n")
                fs.write("\n")
                print("     written to file")
                print("")

    return

def print_avg(avg):
    print(f"  - pressure:           [ {avg['p']:.3e}, 'Pa'     ]")
    print(f"  - temperature:        [ {avg['T']:.3e}, 'K'      ]")
    print(f"  - density:            [ {avg['rho']:.3e}, 'kg/m^3' ]")
    print(f"  - Mach-number:        [ {avg['Mach']:.3e}, '-'      ]")
    print(f"  - x-velocity:         [ {avg['xvel']:.3e}, 'm/s'    ]")
    print(f"  - total-pressure:     [ {avg['p0']:.3e}, 'Pa'     ]")
    print(f"  - total-temperature:  [ {avg['T0']:.3e}, 'K'      ]")
    print(f"  - total-density:      [ {avg['rho0']:.3e}, 'kg/m^3' ]")
    return

def set_gamma_and_r(slice):
    global GAMMA, R
    p = slice.cell_data['p'][0]
    rho = slice.cell_data['rho'][0]
    e = slice.cell_data['e'][0]
    T = slice.cell_data['T'][0]
    GAMMA = p/(rho*e) + 1.0
    Cv = e/T
    R = (GAMMA - 1.0)*Cv
    print(f"# GAMMA set: \t {GAMMA:.2f}")
    print(f"# R set: \t {R:.2f}")
    return

# stagnation expressions taken from gdtk.ideal_gas_flow

def T0_T(M, g=1.4):
    """
    Total to static temperature ratio for an adiabatic flow.

    M: Mach number
    g: ratio of specific heats
    Returns: T0/T
    """
    return 1.0 + (g - 1.0) * 0.5 * M**2

def p0_p(M, g=1.4):
    """
    Total to static pressure ratio for an isentropic flow.

    M: Mach number
    g: ratio of specific heats
    Returns: p0/p
    """
    return (T0_T(M, g))**( g / (g - 1.0) )

def r0_r(M, g=1.4):
    """
    Stagnation to free-stream density ratio for an isentropic flow.

    M: Mach number
    g: ratio of specific heats
    Returns: r0/r
    """
    return (T0_T(M, g))**(1.0 / (g - 1.0))
        
def compute_slice_averages(slice, xLoc):
    area_avg, area = area_weighted_average(slice)
    mass_flow_avg = mass_flow_weighted_average(slice)
    cmes_avg = cmes_average(slice, area, xLoc)
    return area_avg, mass_flow_avg, cmes_avg, area
    

def area_weighted_average(slice):
    avg = dict.fromkeys(Fields, 0.0)
    area = 0.0
    for i in range(slice.n_cells):
        dA = slice.cell_data['area-zy'][i]
        area += dA
        for k in avg.keys():
            avg[k] += slice.cell_data[k][i]*dA
    for k in avg.keys():
        avg[k] /= area
        avg[k] = float(avg[k])
    avg['p0'] = avg['p']*p0_p(avg['Mach'], GAMMA)
    avg['T0'] = avg['T']*T0_T(avg['Mach'], GAMMA)
    avg['rho0'] = avg['rho']*r0_r(avg['Mach'], GAMMA)
    return avg, area

def mass_flow_weighted_average(slice):
    avg = dict.fromkeys(Fields, 0.0)
    f_mass = 0.0
    for i in range(slice.n_cells):
        dA = slice.cell_data['area-zy'][i]
        rho = slice.cell_data['rho'][i]
        u = slice.cell_data['xvel'][i]
        df = rho*u*dA
        f_mass += df
        for k in avg.keys():
            avg[k] += slice.cell_data[k][i]*df
    for k in avg.keys():
        avg[k] /= f_mass
        avg[k] = float(avg[k])
    avg['p0'] = avg['p']*p0_p(avg['Mach'], GAMMA)
    avg['T0'] = avg['T']*T0_T(avg['Mach'], GAMMA)
    avg['rho0'] = avg['rho']*r0_r(avg['Mach'], GAMMA)
    return avg

def compute_fluxes(slice):
    f_mass = 0.0
    f_energy = 0.0
    f_entropy = 0.0
    for i in range(slice.n_cells):
        dA = slice.cell_data['area-zy'][i]
        rho = slice.cell_data['rho'][i]
        u = slice.cell_data['xvel'][i]
        v = slice.cell_data['yvel'][i]
        w = slice.cell_data['zvel'][i]
        p = slice.cell_data['p'][i]
        e = slice.cell_data['e'][i]
        T = slice.cell_data['T'][i]
        h = e + p/rho
        ke = 0.5*(u*u + v*v + w*w)
        h0 = h + ke
        s = ((GAMMA*R)/(GAMMA-1.0))*log(T/T_ref) - R*log(p/p_ref)
        f_mass += (rho*u)*dA
        f_energy += (rho*u*h0)*dA
        f_entropy += (rho*u*s)*dA
    return (f_mass, f_energy, f_entropy)

def cmes_average(slice, area, xLoc):
    f_mass, f_energy, f_entropy = compute_fluxes(slice)
    h0 = f_energy/f_mass
    def fzero(T):
        h = R*T*(1./(GAMMA-1.) + 1.)
        u = sqrt(2.*(h0 - h))
        p = f_mass*R*T/(u*area)
        F = (f_entropy/f_mass) - (((GAMMA*R)/(GAMMA-1.))*log(T/T_ref) - R*log(p/p_ref))
        return F

    x0 = median(slice.cell_data['T'])
    x1 = mean(slice.cell_data['T'])
    sol = root_scalar(fzero, x0=x0, x1=x1, xtol=1.0e-3, rtol=1.0e-6)

    if not sol.converged:
        print(f"pvn2oned: Difficulty finding iterative solution for temperature at x= {xLoc}.")
        print(f"pvn2oned: median T= {x0}, mean T = {x1}, sonic T= {T_sonic}")
        print("pvn2oned: Exiting.")
        sys.exit(1)

    T = sol.root
    # Check T < T_sonic
    h = R*T*(1./(GAMMA-1.) + 1.)
    T_sonic = 2.*(h0 - h)/(GAMMA*R)
    if T > T_sonic:
        print(f"pvn2oned: Iterative solve found temperature value above sonic temperature at x = {xLoc}.")
        print("pvn2oned: This indicates a (predominantly) subsonic flow field.")
        print("pvn2oned: Exiting.")
        sys.exit(1)
    
    # now compute averaged properties
    u = sqrt(2.*(h0 - h))
    p = f_mass*R*T/(u*area)
    rho = p/(R*T)
    a = sqrt(GAMMA*R*T)
    avg = dict.fromkeys(Fields, 0.0)
    avg['rho'] = rho
    avg['p'] = p
    avg['T'] = T
    avg['xvel'] = u
    avg['Mach'] = u/a
    avg['p0'] = avg['p']*p0_p(avg['Mach'], GAMMA)
    avg['T0'] = avg['T']*T0_T(avg['Mach'], GAMMA)
    avg['rho0'] = avg['rho']*r0_r(avg['Mach'], GAMMA)
    return avg

if __name__ == '__main__':
    run()
    
    





