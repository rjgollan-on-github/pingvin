package pingvin

import "core:os"
import "core:fmt"
import "core:math"

write_flow_field_as_vtk :: proc (filename: string) {
    f, err := os.open(filename, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o644)
    defer os.close(f)

    fmt.fprintln(f, "# vtk DataFile Version 2.0")
    fmt.fprintfln(f, "grid_3d, created by pingvin %s", PINGVIN_VERSION)
    fmt.fprintln(f, "ASCII")
    fmt.fprintln(f, "DATASET UNSTRUCTURED_GRID")
    fmt.fprintfln(f, "POINTS %d double", len(global_data.vertices))
    for v in global_data.vertices {
        fmt.fprintfln(f, "%.16e %.16e %.16e", uoflowz(real(v.x)), uoflowz(real(v.y)), uoflowz(real(v.z)))
    }
    n_cells := len(global_data.hexes)
    fmt.fprintln(f, "")
    fmt.fprintfln(f, "CELLS %d %d", n_cells, n_cells*(VtxCount[.hex]+1))
    for h in global_data.hexes {
        fmt.fprintfln(f, "%d %d %d %d %d %d %d %d %d", VtxCount[.hex],
                      h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7])
    }
    fmt.fprintfln(f, "CELL_TYPES %d", n_cells)
    for i in 0..<n_cells {
        fmt.fprintfln(f, "%d", VTKElement.hex)
    }

    fmt.fprintfln(f, "CELL_DATA %d", n_cells)
    for q in Primitive_Quantities {
        fmt.fprintfln(f, "SCALARS %s double", q)
        fmt.fprintln(f, "LOOKUP_TABLE default")
        for cell in global_data.cells {
            fmt.fprintfln(f, "%.18e", real(cell.pqs[q]))
        }
    }
    // Add some variables of interest
    fmt.fprintln(f, "SCALARS Mach double")
    fmt.fprintfln(f, "LOOKUP_TABLE default")
    for cell in global_data.cells {
        speed := math.sqrt(real(cell.pqs[.xvel])*real(cell.pqs[.xvel]) +
            real(cell.pqs[.yvel])*real(cell.pqs[.yvel]) +
            real(cell.pqs[.zvel])*real(cell.pqs[.zvel]))
        a := math.sqrt(real(globals.gamma)*real(cell.pqs[.p])/real(cell.pqs[.rho]))
        M := speed / a
        fmt.fprintfln(f, "%.18e", M)
    }
}

