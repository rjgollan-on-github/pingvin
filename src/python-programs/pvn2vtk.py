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
    # Read in number of points
    nPoints = int(lines[1].split()[1])
    
    vtuDataset = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points_start = 2
    print("Gathering points.")
    for id, line in enumerate(lines[points_start:points_start+nPoints]):
        vtx = list(map(float, line.split()))
        vtx.pop() # drop the index at end
        points.InsertPoint(id, vtx)
    print("Done.")
    print("Setting in dataset")
    vtuDataset.SetPoints(points)
    print("Done.")

    # Get number of cells
    cells_start = points_start+nPoints+1 # actual cells entries
                                         # header is one line before
    nCells = int(lines[cells_start-1].split()[1])
    print("Allocating cells")
    vtuDataset.Allocate(nCells)
    print("Done.")
    print("Gathering cells.")
    for line in lines[cells_start:cells_start+nCells]:
        ids = list(map(int, line.split()))
        ids.pop() # drop id on end
        vtuDataset.InsertNextCell(vtk.VTK_HEXAHEDRON, 8, ids)
    print("Done.")

    del lines


    with gzip.open(FIELD) as f:
        fileContent = f.read().decode()

    lines = fileContent.split('\n')
    vars = lines[0].split()

    for i, var in enumerate(vars):
        print(f"Gathering data for {var}")
        array = vtk.vtkDoubleArray()
        array.SetNumberOfComponents(1)
        array.SetNumberOfTuples(nCells)
        array.SetName(var)
        line_start = 1 + i*nCells
        for id, line in enumerate(lines[line_start:line_start+nCells]):
            val = float(line),
            array.SetTuple(id, val)
        print("   Adding array to dataset")
        vtuDataset.GetCellData().AddArray(array)
        print("   Done.")
        print(f"Done for  {var}.")

    print("Writing...")
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(OUTPUT)
    print("Attaching input data")
    writer.SetInputData(vtuDataset)
    print("Done.")
    print("Writing to disk.")
    writer.Write()
    print("Done.")

    
        


write_vtu_file()

    
    
    

