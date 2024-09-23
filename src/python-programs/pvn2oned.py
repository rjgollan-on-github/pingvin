# Author: RJG
# Date: 2024-09-22

import click
import pyvista as pv

VTUFILE = "pvnsim/flowfield.vtu"

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
        area_avg, mass_flow_avg, area = compute_slice_averages(slice)

        
        print(f"# Averaged properties at x = {xLoc:.3f}")
        print("---")
        print(f"x-location:      [{xLoc:.6e}, 'm']")
        print(f"area-slice:      [{area:.6e}, 'm^2']")
        print(f"number-of-cells:  {slice.n_cells}")
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
            
        
def compute_slice_averages(slice):
    area_avg, area = area_weighted_average(slice)
    mass_flow_avg = mass_flow_weighted_average(slice)
    return area_avg, mass_flow_avg, area
    

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
        



if __name__ == '__main__':
    run()
    
    





