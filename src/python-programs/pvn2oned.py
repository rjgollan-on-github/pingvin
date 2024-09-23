# Author: RJG
# Date: 2024-09-22

import click
import pyvista as pv
from numpy import median
import sys
from math import sqrt
from scipy.optimize import minimize

VTUFILE = "pvnsim/flowfield.vtu"
GAMMA = 1.4
R = 287.1

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
@click.option("-o", "--output-file", "outFile",
              default='pingvin-oned',
              help="base filename for output of slice data",
              show_default=True)
def run(xLoc, nSlices, outFile):

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
        reader = pv.get_reader(VTUFILE)
        mesh = reader.read()
    else:
        print("pvn2oned: No action requested.")
        print("pvn2oned: Consider using one or both of:")
        print("   --at-x-loc X.Y")
        print("   --n-slices N")
        print("pvnoned: Exiting.")
        return

    if printSingleLoc:
        slice = mesh.slice(origin=(xLoc, 0.0, 0.0), normal=(1.0, 0.0, 0.0))
        area_avg, mass_flow_avg, st_avg, area = compute_slice_averages(slice, xLoc)

        
        print(f"# Averaged properties at x = {xLoc:.3f}")
        print("---")
        print(f"x-location:      [{xLoc:.6e}, 'm']")
        print(f"area-slice:      [{area:.6e}, 'm^2']")
        print(f"number-of-cells:  {slice.n_cells}")
        print("")
        print("stream-thrust-average:")
        print_avg(st_avg)
        print("")
        print("mass-flow-weighted-average:")
        print_avg(mass_flow_avg)
        print("")
        print("area-weighted-average:")
        print_avg(area_avg)
        print("...")

    return

def print_avg(avg):
    print(f"  - pressure:       [{avg['p']:.3e}, 'Pa']")
    print(f"  - temperature:    [{avg['T']:.3e}, 'K']")
    print(f"  - density:        [{avg['rho']:.3e}, 'kg/m^3']")
    print(f"  - Mach-number:    [{avg['Mach']:.3e}, '-']")
    print(f"  - x-velocity:     [{avg['xvel']:.3e}, 'm/s']")
    print(f"  - y-velocity:     [{avg['yvel']:.3e}, 'm/s']")
    print(f"  - z-velocity:     [{avg['zvel']:.3e}, 'm/s']")
    return
            
        
def compute_slice_averages(slice, xLoc):
    area_avg, area = area_weighted_average(slice)
    mass_flow_avg = mass_flow_weighted_average(slice)
    st_avg = stream_thrust_average(slice, area, xLoc)
    return area_avg, mass_flow_avg, st_avg, area
    

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
    return avg

def compute_fluxes(slice):
    f_mass = 0.0
    f_xmom = 0.0
    f_ymom = 0.0
    f_zmom = 0.0
    f_energy = 0.0
    for i in range(slice.n_cells):
        dA = slice.cell_data['area-zy'][i]
        rho = slice.cell_data['rho'][i]
        u = slice.cell_data['xvel'][i]
        v = slice.cell_data['yvel'][i]
        w = slice.cell_data['zvel'][i]
        p = slice.cell_data['p'][i]
        e = slice.cell_data['e'][i]
        ke = 0.5*(u*u + v*v + w*w)
        E = e + ke
        f_mass += rho*u*dA
        f_xmom += (rho*u*u + p)*dA
        f_ymom += (rho*u*v)*dA
        f_zmom += (rho*u*w)*dA
        f_energy += (rho*u*E + p*u)*dA
    return (f_mass, f_xmom, f_ymom, f_zmom, f_energy)

def stream_thrust_average(slice, area, xLoc):
    f_mass, f_xmom, f_ymom, f_zmom, f_energy = compute_fluxes(slice)

    def f(x):
        rho, e, u, v, w = x
        p = rho*(GAMMA - 1.0)*e
        E = e + 0.5*(u*u + v*v + w*w)
        # Compute errors
        fmass_err = abs(f_mass - rho*u*area)/(abs(f_mass))
        fxmom_err = abs(f_xmom - (rho*u*u + p)*area)/(abs(f_xmom))
        fymom_err = abs(f_ymom - (rho*u*v)*area)/(abs(f_ymom))
        fzmom_err = abs(f_zmom - (rho*u*w)*area)/(abs(f_zmom))
        fenergy_err = abs(f_energy - (rho*u*E + p*u)*area)/(abs(f_energy))
        err = fmass_err + fxmom_err + fymom_err + fzmom_err + fenergy_err
        return err

    # set guesses for minimizer
    rho = median(slice.cell_data['rho']); rho_min = min(slice.cell_data['rho']); rho_max = max(slice.cell_data['rho'])
    e = median(slice.cell_data['e']); e_min = min(slice.cell_data['e']); e_max = max(slice.cell_data['e'])
    u = median(slice.cell_data['xvel']); u_min = min(slice.cell_data['xvel']); u_max = max(slice.cell_data['xvel'])
    v = median(slice.cell_data['yvel']); v_min = min(slice.cell_data['yvel']); v_max = max(slice.cell_data['yvel'])
    w = median(slice.cell_data['zvel']); w_min = min(slice.cell_data['zvel']); w_max = max(slice.cell_data['zvel'])
    x0 = [rho, e, u, v, w]
    bounds = [(rho_min, rho_max), (e_min, e_max), (u_min, u_max), (v_min, v_max), (w_min, w_max)]
    result = minimize(f, x0, method='Powell', bounds=bounds)

    if not result.success:
        print(f"pvn2oned: Error when attempting to compute stream-thrust average for slice at x= {xLoc}.")
        print("pvn2oned: Median values of slice are: ")
        print(f"  {rho=}  {e=}  {u=}  {v=}  {w=}")
        print("pvn2oned: Best estimate from minimizer: ")
        print(f"  {result.x[0]}  {result.x[1]}  {result.x[2]}  {result.x[3]}  {result.x[4]}")
        print("pvn2oned: Exiting.")
        sys.exit(1)
    
    avg = dict.fromkeys(Fields, 0.0)
    rho, e, u, v, w, = result.x
    avg['rho'] = rho
    p = rho*(GAMMA - 1.0)*e
    avg['p'] = p
    T = p/(rho*R)
    avg['T'] = T
    avg['xvel'] = u
    avg['yvel'] = v
    avg['zvel'] = w
    speed = sqrt(u*u + v*v + w*w)
    a = sqrt(GAMMA*R*T)
    avg['Mach'] = speed/a
    return avg
    
    
    
    

    
        
        
    


if __name__ == '__main__':
    run()
    
    





