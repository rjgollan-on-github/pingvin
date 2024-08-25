package pingvin

import "core:flags"
import "core:os"
import "core:fmt"

test_grid :: "test-assets/qcirc.su2"
test_job_file :: "test-assets/job.lua"
dx :: 0.2

TestSliceCmd := Command {
    main = test_slice,
    description = `[DEV] Test residual evaluation on a slice`,
    short_desc = "''", // ditto
}

test_slice :: proc (args: []string) -> (result: bool) {
    fmt.println("test-slice: Begin...")
    globals.cfg = read_config_from_lua_file(test_job_file)
    read_su2_2d_file(&global_data.up_grid, test_grid)
    free_all(context.temp_allocator)
    allocate_grid_2d(&global_data.dn_grid, len(global_data.up_grid.vertices), len(global_data.up_grid.quads))
    fmt.println("test-slice: Allocated grids.")

    for i in global_data.quads[0] {
        fmt.printfln("i: %d p: %v", i, global_data.vertices[i])
    }

    
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
    apply_inflow(&global_data.slices[0])
    fmt.println("test-slice: Applied inflow conditions to first slice.")

    eval_slice_residual(&global_data.slices[0])
    fmt.println("test-slice: Evaluated residual.")

    delete_global_data()
    return true 
}

