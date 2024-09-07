package pingvin

import "core:os"
import "core:fmt"
import "core:math"

write_flow_field_as_vtk :: proc (filename: string) {
    f, err := os.open(filename, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o644)
    defer os.close(f)
    n_verts := len(global_data.vertices)
    n_cells := len(global_data.cells)

    fmt.fprintln(f, `<VTKFile type="UnstructuredGrid" byte_order="BigEndian">`)
    fmt.fprintln(f, `<UnstructuredGrid>`)
    fmt.fprintfln(f, `<Piece NumberOfPoints="%d" NumberOfCells="%d">`, len(global_data.vertices), len(global_data.cells))

    fmt.fprintln(f, `<Points>`)
    fmt.fprint(f, `  <DataArray type="Float64" NumberOfComponents="3"`)
    fmt.fprintln(f, ` format="ascii">`)
    for v in global_data.vertices {
        fmt.fprintfln(f, "%.18e %.18e %.18e", real(v.x), real(v.y), real(v.z))
    }
    fmt.fprintln(f, `  </DataArray>`)
    fmt.fprintln(f, `</Points>`)

    fmt.fprintln(f, `<Cells>`)
    fmt.fprint(f, `  <DataArray type="Int32" Name="connectivity"`)
    fmt.fprintln(f, ` format="ascii">`)
    for h in global_data.hexes {
        fmt.fprintfln(f, "%d %d %d %d %d %d %d %d", 
                      h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7])
    }
    fmt.fprintln(f, `  </DataArray>`)

    fmt.fprintln(f, `  <DataArray type="Int32" Name="offsets"`)
    fmt.fprintln(f, ` format="ascii">`)
    conn_offset : int = 0;
    for i in 0..<n_cells {
        conn_offset += VtxCount[.hex]
        fmt.fprintfln(f, "%d", conn_offset)
    }
    fmt.fprintln(f, `  </DataArray>`)

    fmt.fprint(f, `  <DataArray type="UInt8" Name="types"`)
    fmt.fprintln(f, ` format="ascii">`)
    for i in 0..<n_cells do fmt.fprintfln(f, "%d", VTKElement.hex)
    fmt.fprintln(f, `  </DataArray>`)
    fmt.fprintln(f, `</Cells>`)

    fmt.fprintln(f, `<CellData>`)
    for q in Primitive_Quantities {
        fmt.fprintf(f, `  <DataArray Name="%s" type="Float64" NumberOfComponents="1"`, q)
        fmt.fprintln(f, ` format="ascii">`)
        for cell in global_data.cells {
            fmt.fprintfln(f, "%.18e", real(cell.pqs[q]))
        }
        fmt.fprintln(f, `  </DataArray>`)
    }
    // Gather velocities as vector
    fmt.fprint(f, `  <DataArray Name="velocity" type="Float64" NumberOfComponents="3"`)
    fmt.fprintfln(f, ` format="ascii">`)
    for cell in global_data.cells {
        fmt.fprintfln(f, `%.18e %.18e %.18e`, real(cell.pqs[.xvel]), real(cell.pqs[.yvel]), real(cell.pqs[.zvel]))
    }
    fmt.fprintln(f, `</DataArray>`)
    // Add some variables of interest
    fmt.fprint(f, `  <DataArray Name="Mach" type="Float64" NumberOfComponents="1"`)
    fmt.fprintln(f, ` format="ascii">`)
    for cell in global_data.cells {
        speed := math.sqrt(real(cell.pqs[.xvel])*real(cell.pqs[.xvel]) +
            real(cell.pqs[.yvel])*real(cell.pqs[.yvel]) +
            real(cell.pqs[.zvel])*real(cell.pqs[.zvel]))
        a := math.sqrt(real(globals.gamma)*real(cell.pqs[.p])/real(cell.pqs[.rho]))
        M := speed / a
        fmt.fprintfln(f, "%.18e", M)
    }
    fmt.fprintln(f, `  </DataArray>`)
    fmt.fprintln(f, `</CellData>`)

    fmt.fprintln(f, `</Piece>`)
    fmt.fprintln(f, `</UnstructuredGrid>`)
    fmt.fprintln(f, `</VTKFile>`)
    return
}

