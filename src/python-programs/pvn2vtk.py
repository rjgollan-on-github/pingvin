# Author: RJG
# Date: 2024-09-08

import gzip
import vtk

GRID = "pvnsim/grid.su2.gz"
FIELD = "pvnsim/flowfield.gz"
OUTPUT = "pvnsim/flowfield.vtu"

def write_vtu_file():
    # 1. Read in grid.
    with gzip.open(GRID) as f:
        fileContent = f.read().decode()

    lines = fileContent.split('\n')
    # Don't need dimensions, always assume 3
    lines.pop(0)
    # Read in number of points
    nPoints = int(lines[0].split()[1])
    lines.pop(0)
    
    vtuDataset = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    for id in range(nPoints):
        vtx = list(map(float, lines[0].split()))
        vtx.pop() # drop the index at end
        points.InsertPoint(id, vtx)
        lines.pop(0)
    vtuDataset.SetPoints(points)

    # Get number of cells
    nCells = int(lines[0].split()[1])
    lines.pop(0)
    vtuDataset.Allocate(nCells)
    for id in range(nCells):
        ids = list(map(int, lines[0].split()))
        lines.pop(0)
        ids.pop() # drop id on end
        vtuDataset.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, ids)

    del lines

    with gzip.open(FIELD) as f:
        fileContent = f.read().decode()

    lines = fileContent.split('\n')
    vars = lines[0].split()
    lines.pop(0)
    for var in vars:
        array = vtk.vtkDoubleArray()
        array.SetNumberOfComponents(1)
        array.SetNumberOfTuples(nCells)
        array.SetName(var)
        for id in range(nCells):
            val = float(lines[0]),
            lines.pop(0)
            array.SetTuple(id, val)
        vtuDataset.GetCellData().AddArray(array)

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(OUTPUT)
    writer.SetInputData(vtuDataset)
    writer.Write()

    
        


write_vtu_file()

    
    
    

