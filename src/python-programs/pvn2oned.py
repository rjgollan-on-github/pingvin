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
R = 287.0

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

    if doSlicing:
        slices = mesh.slice_along_axis(n=nSlices, axis="x", tolerance=0.001)
        vars = ('p', 'T', 'rho', 'Mach', 'xvel', 'yvel', 'zvel')

        with open(f"{outfile}-area-avg.data", "w") as fa, \
             open(f"{outfile}-mass-flow-avg.data", "w") as fm, \
             open(f"{outfile}-stream-thrust-avg.data", "w") as fs:

            header = "x area p T rho Mach xvel yvel zvel\n"
            fa.write(header)
            fm.write(header)
            fs.write(header)
             
            for i, slice in enumerate(slices):
                if i == 0:
                    set_gamma_and_r(slice)
                x = slice.get_cell(0).center[0]
                print(f"-- slice-{i:03d} @ x={x:.3f} --")
                print("     computing averaged properties")
                area_avg, mass_flow_avg, st_avg, area = compute_slice_averages(slice, xLoc)
                x_area_str = f"{x:.6e} {area:.6e} "
                fa.write(x_area_str)
                fm.write(x_area_str)
                fs.write(x_area_str)
                for v in vars:
                    fa.write(f"{area_avg[v]:.6e} ")
                    fm.write(f"{mass_flow_avg[v]:.6e} ")
                    fs.write(f"{st_avg[v]:.6e} ")
                fa.write("\n")
                fm.write("\n")
                fs.write("\n")
                print("     written to file")
                print("")

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
    # set guesses for minimizer
    rho = median(slice.cell_data['rho']); rho_min = min(slice.cell_data['rho']); rho_max = max(slice.cell_data['rho'])
    T = median(slice.cell_data['T']); T_min = min(slice.cell_data['T']); T_max = max(slice.cell_data['T'])
    u = median(slice.cell_data['xvel']); u_min = min(slice.cell_data['xvel']); u_max = max(slice.cell_data['xvel'])
    v = median(slice.cell_data['yvel']); v_min = min(slice.cell_data['yvel']); v_max = max(slice.cell_data['yvel'])
    w = median(slice.cell_data['zvel']); w_min = min(slice.cell_data['zvel']); w_max = max(slice.cell_data['zvel'])

    def f(x):
        rho_n, T_n, u_n, v_n, w_n = x
        rho = rho_n*rho_max
        T = T_n*T_max
        u = u_n*u_max
        v = v_n*(v_max - v_min + 1.0) + v_min
        w = w_n*(w_max - w_min + 1.0) + w_min
        p = rho*R*T
        e = T*R/(GAMMA - 1.0)
        E = e + 0.5*(u*u + v*v + w*w)
        # Compute errors
        fmass_err = abs(f_mass - rho*u*area)/(abs(f_mass))
        fxmom_err = abs(f_xmom - (rho*u*u + p)*area)/(abs(f_xmom))
        fymom_err = abs(f_ymom - (rho*u*v)*area)/(abs(f_ymom))
        fzmom_err = abs(f_zmom - (rho*u*w)*area)/(abs(f_zmom))
        fenergy_err = abs(f_energy - (rho*u*E + p*u)*area)/(abs(f_energy))
        err = fmass_err + fxmom_err + fymom_err + fzmom_err + fenergy_err
        return err

    x0 = [rho/rho_max, T/T_max, u/u_max, (v-v_min)/(v_max-v_min+1.0), (w-w_min)/(w_max-w_min+1.0)]
    bounds = [(rho_min/rho_max, 1.0), (T_min/T_max, 1.0), (u_min/u_max, 1.0), (0.0, 1.0), (0.0, 1.0)]
    result = minimize(f, x0, method='SLSQP', bounds=bounds, options={'ftol':1.0e-6})

    if not result.success:
        print(f"pvn2oned: Error when attempting to compute stream-thrust average for slice at x= {xLoc}.")
        print("pvn2oned: Median values of slice are: ")
        print(f"  {rho=}  {T=}  {u=}  {v=}  {w=}")
        print("pvn2oned: Best estimate from minimizer: ")
        print(f"  {result.x[0]*rho_max}  {result.x[1]*T_max}  {result.x[2]*u_max}  {result.x[3]*(v_max - v_min + 1.0) + v_min}  {result.x[4]*(w_max - w_min + 1.0) + w_min}")
        print("pvn2oned: Exiting.")
        sys.exit(1)
    
    avg = dict.fromkeys(Fields, 0.0)
    rho_n, T_n, u_n, v_n, w_n = result.x
    rho = rho_n*rho_max
    T = T_n*T_max
    u = u_n*u_max
    v = v_n*(v_max - v_min + 1.0) + v_min
    w = w_n*(w_max - w_min + 1.0) + w_min
    avg['rho'] = rho
    p = rho*R*T           
    avg['p'] = p
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
    
    





