package pingvin

import "core:os"
import "core:fmt"
import "core:math"
import "vendor:zlib"
import "core:strings"

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

write_grid_and_field :: proc () {
    if !os.exists("pvnsim") {
        err := os.make_directory("pvnsim")
        if err != os.ERROR_NONE {
            fmt.println("pvn: unable to create directory: 'pvnsim'")
            fmt.println("pvn: exiting")
            os.exit(1)
        }
    }

    grid_name :: `pvnsim/grid.su2.gz`
    write_grid(grid_name)

    field_name :: `pvnsim/flowfield.gz`
    write_flow_field(field_name)
}

write_grid :: proc (filename : cstring) {
    f := zlib.gzopen64(filename, "w0")
    defer zlib.gzclose(f)

    zlib.gzprintf(f, "NDIME= 3\n")

    fmt.printfln("nverts= %d", len(global_data.vertices))
    zlib.gzprintf(f, "NPOIN= %d\n", len(global_data.vertices))
    for v, i in global_data.vertices {
        zlib.gzprintf(f, "%.16e %.16e %.16e %ld\n", f64(real(v.x)), f64(real(v.y)), f64(real(v.z)), i)
    }

    zlib.gzprintf(f, "NELEM= %d\n", len(global_data.hexes))
    for h, i in global_data.hexes {
        zlib.gzprintf(f, "%d %d %d %d %d %d %d %d %d\n",
            h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], i)
    }
    return
}

write_flow_field :: proc (filename : cstring) {
    f := zlib.gzopen64(filename, "w7")
    defer zlib.gzclose(f)

    for q in Primitive_Quantities {
        q_str := fmt.tprintf("%s", q)
        q_cstr := strings.unsafe_string_to_cstring(q_str)
        zlib.gzprintf(f, "%s ", q_cstr)
    }
    zlib.gzprintf(f, "Mach\n")

    for q in Primitive_Quantities {
        for cell in global_data.cells {
            zlib.gzprintf(f, "%.18e\n", f64(real(cell.pqs[q])))
        }
    }
    for cell in global_data.cells {
        speed := math.sqrt(real(cell.pqs[.xvel])*real(cell.pqs[.xvel]) +
            real(cell.pqs[.yvel])*real(cell.pqs[.yvel]) +
            real(cell.pqs[.zvel])*real(cell.pqs[.zvel]))
        a := math.sqrt(real(globals.gamma)*real(cell.pqs[.p])/real(cell.pqs[.rho]))
        M := speed / a
        zlib.gzprintf(f, "%.18e\n", M)
    }
    return
}
