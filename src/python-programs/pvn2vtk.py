# Author: RJG
# Date: 2024-09-08

import click
import gzip
import vtk
import runpy
import sys
import numpy as np
from tqdm import tqdm

GRID = "pvnsim/grid.su2.gz"
FIELD = "pvnsim/flowfield.gz"
OUTPUT = "pvnsim/flowfield.vtu"

@click.command()
@click.option("-q", "--quantities", "quantStr",
              help="comma separated list of variables provided in reference function",
              default=None,
              show_default=True)
@click.option("-r", "--reference-solution", "py_ref_file",
              help="Python file containing reference solution",
              default=None,
              show_default=True)
@click.option("-f", "--reference-function", "ref_fn",
              help="Name of function for reference solution calculation",
              default="ref_soln",
              show_default=True)
def write_vtu_file(quantStr, py_ref_file, ref_fn):
    # 0. Check on function (if necessary)
    if py_ref_file:
        quants = quantStr.split(",")
        quants = [q.strip() for q in quants]
        ref = runpy.run_path(py_ref_file)
        if not ref_fn in ref:
            print(f"pnv2vtk: ERROR -- could not find function {ref_fn} in file {py_ref_file}.")
            print("pvn2vtk: exiting.")
            sys.exit(1)
        # do a test evaluation
        refFn = ref[ref_fn]
        ret = refFn(0.0, 0.0, 0.0)
        if not isinstance(ret, dict):
            print(f"pvn2vtk: ERROR -- reference function {ref_fn} does not return a dict.")
            print("pvn2vtk: exiting.")
            sys.exit(1)
        for q in quants:
            if not q in ret:
                print(f"pvn2vtk: ERROR -- reference function {ref_fn} does not return expected key {q} in dict.")
                print("pvn2vtk: exiting.")
                sys.exit(1)
    # 1. Read in grid.
    with gzip.open(GRID) as f:
        fileContent = f.read().decode()

    lines = fileContent.split('\n')
    # Don't need dimensions, always assume 3
    # Read in number of points
    nPoints = int(lines[1].split()[1])
    
    vtuDataset = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points_start = 2
    for id, line in tqdm(enumerate(lines[points_start:points_start+nPoints]), desc="Gathering points", unit="points", unit_scale=True):
        vtx = list(map(float, line.split()))
        vtx.pop() # drop the index at end
        points.InsertPoint(id, vtx)
    vtuDataset.SetPoints(points)

    # Get number of cells
    cells_start = points_start+nPoints+1 # actual cells entries
                                         # header is one line before
    nCells = int(lines[cells_start-1].split()[1])
    vtuDataset.Allocate(nCells)
    for line in tqdm(lines[cells_start:cells_start+nCells], desc="Gathering cells  ", unit="cells", unit_scale=True):
        ids = list(map(int, line.split()))
        ids.pop() # drop id on end
        vtuDataset.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, ids)

    del lines


    with gzip.open(FIELD) as f:
        fileContent = f.read().decode()

    lines = fileContent.split('\n')
    vars = lines[0].split()

    for i, var in enumerate(vars):
        array = vtk.vtkDoubleArray()
        array.SetNumberOfComponents(1)
        array.SetNumberOfTuples(nCells)
        array.SetName(var)
        line_start = 1 + i*nCells
        for id, line in enumerate(tqdm(lines[line_start:line_start+nCells], desc=f"   Gathering {var} ", unit="cells", unit_scale=True)):
            val = float(line),
            array.SetTuple(id, val)
        vtuDataset.GetCellData().AddArray(array)

    if py_ref_file:
        q_err = dict.fromkeys(quants, vtk.vtkDoubleArray())
        q_ref = dict.fromkeys(quants, vtk.vtkDoubleArray())
        s_arr = dict.fromkeys(quants, 0.0)
        for q in quants:
            s_arr[q] = vtuDataset.GetCellData().GetArray(q)
        for k, a in q_err.items():
            a.SetNumberOfComponents(1)
            a.SetNumberOfTuples(nCells)
            a.SetName(f"{k}-error")
        for k, a in q_ref.items():
            a.SetNumberOfComponents(1)
            a.SetNumberOfTuples(nCells)
            a.SetName(f"{k}-ref")
        centroid = np.zeros(3)
        for i in tqdm(range(nCells), desc=f"Computing reference quantities ", unit="cells", unit_scale=True):
            cell = vtuDataset.GetCell(i)
            cell.GetCentroid(centroid)
            x, y, z = centroid
            vals = refFn(x, y, z)
            for q in quants:
                q_ref[q].SetTuple(i, (vals[q],))
                s_val = s_arr[q].GetValue(i)
                q_err[q].SetTuple(i, (vals[q] - s_val,))
        for k, a in q_err.items():
            vtuDataset.GetCellData().AddArray(a)
        for k, a in q_ref.items():
            vtuDataset.GetCellData().AddArray(a)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(OUTPUT)
    writer.SetInputData(vtuDataset)
    print("Writing to disk.")
    writer.Write()
    print("Done.")

write_vtu_file()

    
    
    

