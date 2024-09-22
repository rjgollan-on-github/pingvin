# Author: RJG
# Date: 2024-09-22

import click
import pyvista as pv

VTUFILE = "pvnsim/flowfield.vtu"

Properties = ['u', 'M', 'p', 'T', 'rho', 'p_tot', 'T_tot']


@click.command()
@click.option("-@", "--at-x-loc", "xLoc",
              help="x location for extraction of single set of one-d properties",
              show_default=False)
@click.option("-n", "--n-slices", "nSlices",
              default=-1,
              help="number of slices for extracting one-d properties (distributed equally in x)",
              show_default=False)
@click.option("-o", "--output-file", "outFile",
              default='pingvin-oned.data',
              help="name of output file for slice data",
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
        (area_avg, ) = compute_slice_averages(slice)
        print(f"Area-weighted average properties at x = {xLoc:.3e}")
        print(area_avg)

    return
    
        
def compute_slice_averages(slice):
    return area_weighted_average(slice), 
    

def area_weighted_average(slice):
    avg = dict.fromkeys(Properties, 0.0)
    area = 0.0
    for i in range(slice.n_cells):
        dA = slice.cell_data['area-zy'][i]
        area += dA
        avg['u'] += slice.cell_data['xvel'][i]*dA
        avg['M'] += slice.cell_data['Mach'][i]*dA
        avg['p'] += slice.cell_data['p'][i]*dA
        avg['T'] += slice.cell_data['T'][i]*dA
        avg['rho'] += slice.cell_data['rho'][i]*dA
    for k in avg.keys():
        avg[k] /= area
        avg[k] = float(avg[k])
    return avg

if __name__ == '__main__':
    run()
    
    





