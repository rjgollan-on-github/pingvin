# Author: RJG
# Date: 2024-09-08

import gzip
import vtk
from tqdm import tqdm

GRID = "pvnsim/grid.su2.gz"
FIELD = "pvnsim/flowfield.gz"
OUTPUT = "pvnsim/flowfield.vtu"

def write_vtu_file():
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
        for id, line in enumerate(tqdm(lines[line_start:line_start+nCells], desc=f"Gathering {var} ")):
            val = float(line),
            array.SetTuple(id, val)
        vtuDataset.GetCellData().AddArray(array)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(OUTPUT)
    writer.SetInputData(vtuDataset)
    print("Writing to disk.")
    writer.Write()
    print("Done.")

write_vtu_file()

    
    
    

