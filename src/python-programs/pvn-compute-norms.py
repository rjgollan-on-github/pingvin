# Author: RJG
# Date: 2025-02-15
#
#

import click
import pyvista as pv
import sys
import runpy
from math import fabs, sqrt
from contextlib import nullcontext
from tqdm import tqdm

VTUFILE = "pvnsim/flowfield.vtu"

@click.command()
@click.option("-n", "--norms", "normStr",
              help="comma separated list of variables for norms calculation",
              default="Mach",
              show_default=True)
@click.option("-r", "--reference-solution", "py_ref_file",
              help="Python file containing reference solution",
              default="reference-solution.py",
              show_default=True)
@click.option("-f", "--reference-function", "ref_fn",
              help="Name of function for reference solution calculation",
              default="ref_soln",
              show_default=True)
@click.option("-s", "--selection-function", "sel_fn",
              help="name of function for determining extents of norms calculation",
              default=None)
@click.option("-o", "--output-file", "outfile",
              help="name of output file",
              default=None,
              show_default=False)
def run(normStr, py_ref_file, ref_fn, sel_fn, outfile):

    norms = normStr.split(",")
    norms = [norm.strip() for norm in norms]

    L1, L2, Linf, pos, vol_avg, n_cells, n_cells_in_calc = compute_norms(norms, py_ref_file, ref_fn, sel_fn)

    with open(outfile, "w") if outfile else nullcontext(sys.stdout) as f:
        f.write("---\n")
        f.write(f"average-volume:    [{vol_avg:.6e}, 'm^3']\n")
        f.write(f"number-cells:      {n_cells}\n")
        f.write(f"number-cells-in-selection: {n_cells_in_calc}\n")
        for norm in norms:
            f.write(f"norm-{norm}:\n")
            f.write(f"  L1:   {L1[norm]:.6e}\n")
            f.write(f"  L2:   {L2[norm]:.6e}\n")
            f.write(f"  Linf: {Linf[norm]:.6e}\n")
            f.write(f"  pos:  {{x: {pos[0]:.6e}, y: {pos[1]:.6e}, z: {pos[2]:.6e}}}\n")
        f.write("...\n")

    return

def compute_norms(norms, py_ref_file, ref_fn, sel_fn):
    pvn = pv.read(VTUFILE)
    ref = runpy.run_path(py_ref_file)
    L1 = dict.fromkeys(norms, 0.0)
    L2 = dict.fromkeys(norms, 0.0)
    Linf = dict.fromkeys(norms, 0.0)
    ctrs = pvn.cell_centers()
    vol_sum = 0.0
    n_cells_in_calc = 0
    for i in tqdm(range(pvn.n_cells), desc="Computing norms", unit="cells", unit_scale=True):
        x, y, z = ctrs.points[i]
        if sel_fn:
            if not ref[sel_fn](x, y, z):
                continue
        vals = ref[ref_fn](x, y, z)
        vol = pvn.cell_data['vol'][i]
        vol_sum += vol
        n_cells_in_calc += 1
        for var in norms:
            vals[var] = pvn.cell_data[var][i] - vals[var]
        for var in norms:
            L1[var] += vol * fabs(vals[var])
            L2[var] += vol * vals[var]**2
            if (fabs(vals[var]) > Linf[var]):
                Linf[var] = fabs(vals[var])
                pos = (x, y, z)
    for var in norms:
        L1[var] /= vol_sum
        L2[var] = sqrt(L2[var]/vol_sum)

    return L1, L2, Linf, pos, vol_sum/n_cells_in_calc, pvn.n_cells, n_cells_in_calc

if __name__ == '__main__':
    run()
    
