package pingvin

import "core:flags"
import "core:os"
import "core:fmt"
import "core:math/rand"

@(private="file")
test_grid :: "test-assets/single-slice/qcirc.su2"
@(private="file")
test_job_file :: "test-assets/single-slice/job.lua"
@(private="file")
dx :: 0.1

TestSliceCmd := Command {
    main = test_slice,
    description = `[DEV] Test Newton solver on a single slice`,
    short_desc = "''", // ditto
}

Initialisation :: enum {freestream_p, rand}

TestSliceOptions :: struct {
    init : Initialisation `args:"name=init" usage:Select initialisation from 'freestream_p' or 'rand'`,
}

test_slice :: proc (args : []string) -> (result : bool) {

    cmd_name := args[0]
    opt := TestSliceOptions{init=.freestream_p}
    err := flags.parse(&opt, args[1:])

    if err != nil {
        flags.print_errors(typeid_of(TestSliceOptions), err, cmd_name)
        os.exit(1)
    }

    fmt.println("test-slice: Begin...")
    globals.cfg = read_config_from_lua_file(test_job_file)
    read_su2_2d_file(&global_data.up_grid, test_grid, 0.0)
    free_all(context.temp_allocator)
    allocate_grid_2d(&global_data.dn_grid, len(global_data.up_grid.vertices), len(global_data.up_grid.quads))
    fmt.println("test-slice: Allocated grids.")

    nvars := len(global_data.up_grid.quads) * len(Conserved_Quantities)
    allocate_GMRES_Workspace(nvars, globals.cfg.max_gmres_iterations)
    fmt.println("test-slice: Allocated GMRES workspace")

    // Make new vertices, but offset  x+dx
    n_verts := len(global_data.vertices)
    for i in 0..<n_verts {
        v := global_data.vertices[i]
        append(&global_data.vertices, Vector3{v.x+dx, v.y, v.z})
    }
    fmt.println("test-slice: Created new vertices.")
    // Make new quads (both globally and locally)
    n_quads := len(global_data.quads)
    n_offset := VtxId(n_verts)
    for i in 0..<n_quads {
        global_data.dn_grid.quads[i] = global_data.up_grid.quads[i] + n_offset
        append(&global_data.quads, global_data.dn_grid.quads[i])
    }
    
    fmt.println("test-slice: Created new quads at downstream location.")
    // Now make 3D slice
    assemble_initial_upstream_interfaces(global_data.up_grid)
    fmt.println("test-slice: Created upstream interfaces.")
    append(&global_data.slices, Slice{})
    slice := Slice_Id(0)
    assemble_slice_cells_and_interfaces(&global_data.slices[slice], global_data.up_grid.quads[:], global_data.dn_grid.quads[:], slice)
    fmt.println("test-slice: Created slice of cells.")

    switch opt.init {
    case .freestream_p:
        apply_inflow(&global_data.slices[0])
        apply_xvel_increment(&global_data.slices[0])
        fmt.println("test-slice: Applied inflow conditions with 10% increase on x-vel to first slice.")
    case .rand:
        fmt.println("test-slice: Applied random perturbations to y- and z-momenta.")
        apply_inflow(&global_data.slices[0])
        apply_random_perturbations(&global_data.slices[0])
    }
    
    is_converged := solve_slice(0)
    fmt.printfln("test-slice: solver converged? %v", is_converged)

    fmt.println("test-slice: Write flow field to VTK")
    write_flow_field_as_vtk("single-slice.vtu")

    fmt.println("test-slice: Write SU2 grid")
    write_su2_grid("single-slice.su2")
    

    fmt.println("test-slice: Attempting to delete global data")
    delete_global_data()
    fmt.println("test-slice: Attempting to delete GMRES workspace")
    delete_GMRES_Workspace()
    fmt.println("test-slice: memory deleted")
    return true 
}

apply_xvel_increment :: proc (slice : ^Slice) {
    for &cell in global_data.cells[slice.first_cell:slice.last_cell+1] {
        cell.pqs[.xvel] *= 1.1
        cell.cqs = cq_from_prim(cell.pqs)
    }
}

apply_random_perturbations :: proc (slice : ^Slice) {
    rand.reset(17)
    base_vel := 1.0
    for &cell in global_data.cells[slice.first_cell:slice.last_cell+1] {
        cell.pqs[.yvel] = complex(base_vel*rand.float64_range(-1, 1), 0.0)
        cell.pqs[.zvel] = complex(base_vel*rand.float64_range(-1, 1), 0.0)
        cell.cqs = cq_from_prim(cell.pqs)
    }
    return
}

write_su2_grid :: proc (filename : string) {
    f, err := os.open(filename, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o644)
    defer os.close(f)

    fmt.fprintln(f, "NDIME= 3")

    fmt.fprintfln(f, "NPOIN= %d", len(global_data.vertices))
    for v, i in global_data.vertices {
        fmt.fprintfln(f, "%.16e %.16e %.16e %d", v.x, v.y, v.z, i)
    }

    fmt.fprintfln(f, "NELEM= %d", len(global_data.hexes))
    for h, i in global_data.hexes {
        fmt.fprintfln(f, "%d %d %d %d %d %d %d %d %d %d", VTKElement.hex,
            h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], i)
    }

    fmt.fprintln(f, "NMARK= 4")
    fmt.fprintln(f, "MARKER_TAG= wall")
    fmt.fprintfln(f, "MARKER_ELEMS= %d", len(global_data.slices[0].wall_faces))
    for f_id in global_data.slices[0].wall_faces {
        q := global_data.x_faces[f_id].quad
        fmt.fprintfln(f, "%d %d %d %d %d", VTKElement.quad, q[0], q[1], q[2], q[3])
    }

    fmt.fprintln(f, "MARKER_TAG= symm")
    fmt.fprintfln(f, "MARKER_ELEMS= %d", len(global_data.slices[0].symm_faces))
    for f_id in global_data.slices[0].symm_faces {
        q := global_data.x_faces[f_id].quad
        fmt.fprintfln(f, "%d %d %d %d %d", VTKElement.quad, q[0], q[1], q[2], q[3])
    }

    fmt.fprintln(f, "MARKER_TAG= inflow")
    fmt.fprintfln(f, "MARKER_ELEMS= %d", len(global_data.slices[0].up_faces))
    for face in global_data.slices[0].up_faces {
        q := face.quad
        fmt.fprintfln(f, "%d %d %d %d %d", VTKElement.quad, q[0], q[1], q[2], q[3])
    }
    
    fmt.fprintln(f, "MARKER_TAG= outflow")
    fmt.fprintfln(f, "MARKER_ELEMS= %d", len(global_data.slices[0].dn_faces))
    for face in global_data.slices[0].dn_faces {
        q := face.quad
        fmt.fprintfln(f, "%d %d %d %d %d", VTKElement.quad, q[0], q[1], q[2], q[3])
    }
    
    return
}

