package pingvin


import "core:os"
import "core:strings"
import "core:fmt"
import "core:strconv"
import "core:flags"
import "core:log"
import "core:math"


SMALL_THETA_TOL := 1.0e-6

Grid_parameterisation :: enum {rtheta, bbox}

/*
 * Errors related to grids
 */
GridInitError :: union {
    FileOpenFailed,
    GridReadingError,
    GridFormationError,
}

// Two boundary types are possible: wall and symm.
// The wall represencts the surface of interest, eg. an inlet or wall of a nozzle.
// 'symm' stands for symmetry. It is a numerical boundary condition to exploit
// symmetry in the domain of interest.
Boundary :: enum{wall, symm}

// A 2D boundary is represented by the collection of faces that make up the boundary.
// Since we are working unstructured, these faces may not join contiguously.
Boundary_2d :: struct {
    faces: [dynamic]Face2,
}

/*
 * Grid types:
 *   Grid_2d
 *   Grid_rtheta
 *   Grid_bbox
 */

// A 2D grid is used to represent planes.
// We work in the z-y plane for 2D grids (since we march in x-direction as the third dimension).
// We need to read in one of these 2D grids at initialisation.
// We build others as we march through the domain.
// New 2d grids are created based on an r-theta base grid.
Grid_2d :: struct {
    vertices: [dynamic]VtxId,
    quads: [dynamic]Quad,
    wall_boundary: ^Boundary_2d,
    symm_boundary: ^Boundary_2d,
}

// A grid defined by radius and angle (r-theta).
// The r value is normalised based on the local radius to the boundary.
// This is used as a base grid for transformation from r-theta space
// to z-y plane positions.
Grid_rtheta :: struct {
    r_bar: []f64,
    theta: []f64,
}

// A grid defined by parameters (u,v) that vary between
// 0 --> 1 with (1,1) located at the upper-right of
// a bounding box for the grid.
Grid_bbox :: struct {
    u: []f64,
    v: []f64,
}


/*
 * Grid management functions
 */
allocate_grid_2d :: proc (grid: ^Grid_2d, n_verts, n_quads: int) {
    grid.vertices = make([dynamic]VtxId, n_verts)
    grid.quads = make([dynamic]Quad, n_quads)
}

delete_grid_2d :: proc(grid: ^Grid_2d) {
    delete(grid.vertices)
    delete(grid.quads)
}

copy_grid_2d :: proc(dst, src: ^Grid_2d) {
    copy(dst.vertices[:], src.vertices[:])
    copy(dst.quads[:], src.quads[:])
}

/*
 * SU2 reading functions, and related helper functions.
 */

/*
Finds markers of the form NTHING= X in SU2 files.

Inputs:
- lines: the lines (split) from the SU2 file
- marker: the string marker to seek

Returns:
- line_start: the *next* line after the marker (because typically we need to pick up processing from that point)
- result: A boolean indicating success or failure at finding 'marker'
*/
find_marker :: proc (lines: []string, marker: string) -> (line_start, n_things: int, result: bool) {
    for line, line_no in lines {
        if strings.has_prefix(line, marker) {
            line_start = line_no + 1
            tokens := strings.fields(line, context.temp_allocator)
            n_things = strconv.atoi(tokens[1])
            return line_start, n_things, true
        }
    }
    // If we get here, error
    return line_start, n_things, false
}

/*
Finds marker tags of the form MARKER_TAG= <somstring> in SU2 files.

Inputs:
- lines: the lines (split) from the SU2 file
- tag: the string tag to seek

Returns:
- line_start: the line number *two after* after the marker_tag (because typically we need to pick up processing from that point)
- n_elems: how many elements associated with the marker tag
- result: A boolean indicating success or failure at finding 'marker_tag'
*/
find_marker_tag :: proc (lines: []string, tag: string) -> (line_start, n_elems: int, result: bool) {
    match := fmt.tprintf("MARKER_TAG= %s", tag)
    for line, line_no in lines {
        if strings.compare(line, match) == 0 {
            line_start = line_no + 2
            tokens := strings.fields(lines[line_no+1], context.temp_allocator)
            n_elems = strconv.atoi(tokens[1])
            return line_start, n_elems, true
        }
    }
    return line_start, n_elems, false
}

read_su2_2d_file :: proc (grid: ^Grid_2d, filepath: string) -> (err: GridInitError) {
    data, ok := os.read_entire_file(filepath)
    if !ok {
        return FileOpenFailed{filename=filepath}
    }
    defer delete(data)

    lines := strings.split_lines(string(data), context.temp_allocator)
    // 0. Examine dimensions on first line
    tokens := strings.fields(lines[0], context.temp_allocator)
    dim := strconv.atoi(tokens[1])
    if dim != 2 {
        msg := fmt.tprintf("Error in SU2 file '%s'. NDIME= 2 expected.", filepath);
        return GridFormationError{msg=msg}
    }
    
    // 1. Find points
    line_points_start, n_points : int
    line_points_start, n_points, ok = find_marker(lines, "NPOIN")
    if !ok {
        msg := fmt.tprintf("Error reading NPOIN in SU2 file '%s'", filepath)
        return GridReadingError{msg=msg}
    }
    // Initialise vertices array and start to populate
    global_data.vertices = make([dynamic]Vector3, n_points)
    for line in lines[line_points_start:][:n_points] {
        tokens = strings.fields(line, context.temp_allocator)
        if len(tokens) != 3 {
            msg := fmt.tprintf("Error reading points in SU2 file '%s'. Three items expected on line.", filepath)
            return GridReadingError{msg=msg}
        }
        z := strconv.atof(tokens[0])
        y := strconv.atof(tokens[1])
        i := strconv.atoi(tokens[2])
        global_data.vertices[i] = {complex(0.0, 0.0), complex(y, 0.0), complex(z, 0.0)}
        append(&(grid.vertices), VtxId(i))
    }

    // 2. Find elements
    line_elem_start, n_elem : int
    line_elem_start, n_elem, ok = find_marker(lines, "NELEM")
    if !ok {
        msg := fmt.tprintf("Error reading NELEM in SU2 file '%s'", filepath)
        return GridReadingError{msg=msg}
    }
    // This looks like duplication, but...
    // We want to store a single global list of
    // quad integers because this is useful to
    // write out 2D slices if needed (without recomputing integers
    // for a VTK files).
    // However, we also need grid-specific access of quads to the
    // global vertices array.
    global_data.quads = make([dynamic]Quad, n_elem)
    grid.quads = make([dynamic]Quad, n_elem)
    for line in lines[line_elem_start:][:n_elem] {
        tokens = strings.fields(line, context.temp_allocator)
        if len(tokens) != 6 {
            msg := fmt.tprintf("Error reading quad in SU2 file '%s'. Six items expected on line.", filepath)
            return GridReadingError{msg=msg}
        }
        e_type := strconv.atoi(tokens[0])
        if (e_type != int(VTKElement.quad)) {
            msg := fmt.tprintf("Error in SU2 file '%s'. Quad element expected (9); found %d.", filepath, e_type)
            return GridReadingError{msg=msg}
        }
        i := strconv.atoi(tokens[5])
        // Assemble quad in what appears to be reverse order.
        // The reason is that the GMSH grid was built in the x-y plane
        // We now want to work in z-y plane, so this reversal gives
        // us a face normal that points in positive x direction.
        global_data.quads[i] = {VtxId(strconv.atoi(tokens[4])),
                                VtxId(strconv.atoi(tokens[3])),
                                VtxId(strconv.atoi(tokens[2])),
                                VtxId(strconv.atoi(tokens[1]))}
        grid.quads[i] = global_data.quads[i]
    }

    // 3. Collate boundaries
    line_mark_start, n_mark : int
    _, n_mark, ok = find_marker(lines, "NMARK")
    if !ok {
        msg := fmt.tprintf("Error reading NMARK in SU2 file '%s'", filepath)
        return GridReadingError{msg=msg}
    }

    if n_mark != 2 {
        msg := fmt.tprintf("Error in SU2 file. Exactly two marker tags expected. Found %d.", n_mark)
        return GridFormationError{msg=msg}
    }
    
    // 3a. First grab wall.
    line_wall_elems, n_wall_elems : int
    line_wall_elems, n_wall_elems, ok = find_marker_tag(lines, "wall")
    if !ok {
        msg := fmt.tprintf("Error in SU2 file. 'MARKER_TAG= wall' not found.")
        return GridFormationError{msg=msg}
    }
    global_data.wall_boundary.faces = make([dynamic]Face2, n_wall_elems)
    grid.wall_boundary = &(global_data.wall_boundary)
    for line, i in lines[line_wall_elems:][:n_wall_elems] {
        tokens = strings.fields(line, context.temp_allocator)
        if len(tokens) != 3 {
            msg := fmt.tprintf("Error reading quad in SU2 file '%s'. Three items expected on line.", filepath)
            return GridReadingError{msg=msg}
        }
        e_type := strconv.atoi(tokens[0])
        if (e_type != int(VTKElement.line)) {
            msg := fmt.tprintf("Error in SU2 file '%s'. Line element expected (3); found %d.", filepath, e_type)
            return GridReadingError{msg=msg}
        }
        grid.wall_boundary.faces[i] = {VtxId(strconv.atoi(tokens[1])), VtxId(strconv.atoi(tokens[2]))}
    }
    
    // 3b. Next grab symm.
    line_symm_elems, n_symm_elems : int
    line_symm_elems, n_symm_elems, ok = find_marker_tag(lines, "symm")
    if !ok {
        msg := fmt.tprintf("Error in SU2 file. 'MARKER_TAG= symm' not found.")
        return GridFormationError{msg=msg}
    }
    global_data.symm_boundary.faces = make([dynamic]Face2, n_symm_elems)
    grid.symm_boundary = &(global_data.symm_boundary)
    for line, i in lines[line_symm_elems:][:n_symm_elems] {
        tokens = strings.fields(line, context.temp_allocator)
        if len(tokens) != 3 {
            msg := fmt.tprintf("Error reading quad in SU2 file '%s'. Three items expected on line.", filepath)
            return GridReadingError{msg=msg}
        }
        e_type := strconv.atoi(tokens[0])
        if (e_type != int(VTKElement.line)) {
            msg := fmt.tprintf("Error in SU2 file '%s'. Line element expected (3); found %d.", filepath, e_type)
            return GridReadingError{msg=msg}
        }
        grid.symm_boundary.faces[i] = {VtxId(strconv.atoi(tokens[1])), VtxId(strconv.atoi(tokens[2]))}
    }

    return nil
}

/*
 * Writing functions related to writing grids to VTK format.
 */

write_grid_2d_as_vtk :: proc (filename: string, g: ^Grid_2d) {
    f, err := os.open(filename, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o644)
    defer os.close(f)

    fmt.fprintln(f, "# vtk DataFile Version 2.0")
    fmt.fprintfln(f, "grid_2d, created by pingvin %s", PINGVIN_VERSION)
    fmt.fprintln(f, "ASCII")
    fmt.fprintln(f, "DATASET UNSTRUCTURED_GRID")
    fmt.fprintfln(f, "POINTS %d double", len(g.vertices))
    for v_id in g.vertices {
        v := global_data.vertices[v_id]
        fmt.fprintfln(f, "%.16e %.16e %.16e", real(v.x), real(v.y), real(v.z))
    }
    n_cells := len(global_data.quads)
    fmt.fprintln(f, "")
    fmt.fprintfln(f, "CELLS %d %d", n_cells, n_cells*(VtxCount[.quad]+1))
        for q in global_data.quads {
        fmt.fprintfln(f, "%d %d %d %d %d", VtxCount[.quad], q[0], q[1], q[2], q[3])
    }
    fmt.fprintfln(f, "CELL_TYPES %d", n_cells)
    for i in 0..<n_cells {
        fmt.fprintfln(f, "%d", VTKElement.quad)
    }
    
}

magnitude_yz :: proc (v: Vector3) -> f64 {
    return math.sqrt(real(v.z)*real(v.z) + real(v.y)*real(v.y))
}

argument_yz :: proc (v: Vector3) -> f64 {
    return math.atan2(real(v.y), real(v.z))
}

find_cross_section_points :: proc (xsect: ^Cross_Section, theta: f64) -> (pA, pB: Vector3, result: bool) {
    for i in 0..<len(xsect.vertices)-1 {
        p0 := xsect.vertices[i]
        p1 := xsect.vertices[i+1]
        theta0 := math.atan2(real(p0.y), real(p0.z))
        theta1 := math.atan2(real(p1.y), real(p1.z))
        if (theta0 <= theta) && (theta <= theta1) {
            pA = p0
            pB = p1
            return pA, pB, true
        }
    }
    return pA, pB, false
}

compute_rtheta_grid :: proc (rtg: ^Grid_rtheta, g: ^Grid_2d, xsect: ^Cross_Section) {
    for vtxId, i in g.vertices {
        v := global_data.vertices[vtxId]
        r := magnitude_yz(v)
        theta := argument_yz(v)
        pA, pB, ok := find_cross_section_points(xsect, theta)
        if !ok {
            fmt.printfln("Error in compute_rtheta_grid: r= %v theta= %v", r, theta)
            os.exit(1)
        }
        theta0 := argument_yz(pA)
        theta1 := argument_yz(pB)
        w := (theta - theta0)/(theta1 - theta0)
        r_b := (1 - w)*magnitude_yz(pA) + w*magnitude_yz(pB)
        // Now populate the rtheta grid
        rtg.r_bar[i] = r/r_b
        rtg.theta[i] = theta
    }
}

compute_bbox_grid :: proc (bbg: ^Grid_bbox, g: ^Grid_2d, corner: Vector3) {
    z_c := real(corner.z)
    y_c := real(corner.y)
    for vtxId, i in g.vertices {
        v := global_data.vertices[vtxId]
        z := real(v.z)
        y := real(v.y)

        bbg.u[i] = z/z_c
        bbg.v[i] = y/y_c
    }
}

compute_grid_2d :: proc (g, g_prev: ^Grid_2d, rtg: ^Grid_rtheta, xsect: ^Cross_Section) {
    x := real(xsect.vertices[0].x)
    for i in 0..<len(rtg.r_bar) {
        theta := rtg.theta[i]
        //fmt.printfln("looking for theta= %.6f", theta)
        pA, pB, ok := find_cross_section_points(xsect, theta)
        w := -1.0 // Set negative as sentinel
        if !ok {
            // Catch some end point cases
            p0 := xsect.vertices[0]
            theta0 := math.atan2(real(p0.y), real(p0.z))
            if (theta0 - theta) < SMALL_THETA_TOL {
                pA = p0
                w = 0.0
            }
            p1 := xsect.vertices[len(xsect.vertices)-1]
            theta1 := math.atan2(real(p1.y), real(p1.z))
            if (theta - theta1) < SMALL_THETA_TOL {
                pA = p1
                w = 0.0
            }

            if w != 0.0 {
                fmt.printfln("Error in compute_grid_2d: theta= %v", theta)
                os.exit(1)
            }
        }
        r_b, r, z, y : f64
        if w != 0.0 {
            theta0 := math.atan2(real(pA.y), real(pA.z))
            theta1 := math.atan2(real(pB.y), real(pB.z))
            w = (theta - theta0)/(theta1 - theta0)
            r_b = (1 - w)*magnitude_yz(pA) + w*magnitude_yz(pB)
        }
        else {
            r_b = magnitude_yz(pA)
        }
        r = rtg.r_bar[i]*r_b
        z = r*math.cos(theta)
        y = r*math.sin(theta)
        append(&(global_data.vertices), Vector3{complex(x, 0.0), complex(y, 0.0), complex(z, 0.0)})
        g.vertices[i] = VtxId(len(global_data.vertices)-1)
    }
    n_offset := VtxId(len(g.vertices))
    for i in 0..<len(g.quads) {
        g.quads[i] = g_prev.quads[i] + n_offset
    } 
}

compute_grid_2d_from_bbox :: proc (g, g_prev: ^Grid_2d, bbg: ^Grid_bbox, corner: Vector3) {
    x := real(corner.x)
    for i in 0..<len(bbg.u) {
        u := bbg.u[i]
        v := bbg.v[i]

        z := u*real(corner.z)
        y := v*real(corner.y)

        append(&(global_data.vertices), Vector3{complex(x, 0.0), complex(y, 0.0), complex(z, 0.0)})
        g.vertices[i] = VtxId(len(global_data.vertices)-1)
    }
    n_offset := VtxId(len(g.vertices))
    for i in 0..<len(g.quads) {
        g.quads[i] = g_prev.quads[i] + n_offset
    }
}

/*
 * Functions for building 3D cells from 2D grid planes.
 */

hex_from_quads :: proc (q0, q1: Quad) -> (hex: Hex) {
    hex = Hex{q0[0], q0[1], q0[2], q0[3], 
              q1[0], q1[1], q1[2], q1[3]}
    return hex
}

add_3d_slice_of_hexes :: proc (gu, gd: ^Grid_2d) -> (result: bool) {
    for i in 0..<len(gu.quads) {
        hex := hex_from_quads(gu.quads[i], gd.quads[i])
        append(&global_data.hexes, hex)
    }
    return true
}

/*
 * Functions for writing 3D representation of grid
 */

write_grid_3d_as_vtk :: proc (filename: string) {
    f, err := os.open(filename, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o644)
    defer os.close(f)

    fmt.fprintln(f, "# vtk DataFile Version 2.0")
    fmt.fprintfln(f, "grid_3d, created by pingvin %s", PINGVIN_VERSION)
    fmt.fprintln(f, "ASCII")
    fmt.fprintln(f, "DATASET UNSTRUCTURED_GRID")
    fmt.fprintfln(f, "POINTS %d double", len(global_data.vertices))
    for v in global_data.vertices {
        fmt.fprintfln(f, "%.16e %.16e %.16e", real(v.x), real(v.y), real(v.z))
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
    
}

/*
 * Functions for building the complete 3D grid.
 */

generate_3d_grid :: proc (cfg: Config) -> (result: bool) {
    start : f64
    end : f64
    switch cfg.grid_parameterisation {
    case .rtheta:
        start = real(global_data.xsects[0].vertices[0].x)
        end = real(global_data.xsects[len(global_data.xsects)-1].vertices[0].x)
    case .bbox:
        start = real(global_data.bbox.corners[0].x)
        end = real(global_data.bbox.corners[len(global_data.bbox.corners)-1].x)
    }
    dx := cfg.dx

    // Read grid at first plane
    up_grid : Grid_2d
    err := read_su2_2d_file(&up_grid, cfg.grid2d_file)
    if err != nil {
        fmt.printfln("pvn: Error reading 2D grid file '%s'", cfg.grid2d_file)
    }
    dn_grid : Grid_2d
    allocate_grid_2d(&dn_grid, len(up_grid.vertices), len(up_grid.quads))
    defer delete_grid_2d(&up_grid)
    defer delete_grid_2d(&dn_grid)

    loft_end : f64
    idx_loft_end : int
    curr_loft : Cross_Section_Loft
    curr_xsect : Cross_Section
    curr_rail : Bbox_rail

    switch cfg.grid_parameterisation {
    case .rtheta:
        // Prepare (global) r-theta grid
        allocate_rtheta_grid(len(up_grid.vertices))
        compute_rtheta_grid(&global_data.rtheta_grid, &up_grid, &global_data.xsects[0])

        // Create initial loft section
        n_seg := len(global_data.xsects[0].vertices)
        allocate_cross_section_loft(&curr_loft, n_seg)

        create_cross_section_loft(&curr_loft, &global_data.xsects[0], &global_data.xsects[1])
        loft_end = real(global_data.xsects[1].vertices[0].x)
        idx_loft_end = 1

        allocate_cross_section(&curr_xsect, n_seg)
    case .bbox:
        allocate_bbox_grid(len(up_grid.vertices))
        compute_bbox_grid(&global_data.bbox_grid, &up_grid, global_data.bbox.corners[0])

        // Create initial rail
        curr_rail = create_bbox_rail(&global_data.bbox, 1)
        loft_end = real(curr_rail.end.x)
        idx_loft_end = 1
    }

    for x := start + dx; x < (end + 0.1*dx); x += dx {
        // We might need to create a new loft
        if x > loft_end {
            switch cfg.grid_parameterisation {
            case .rtheta:
                // Search for new loft end in cross sections, beginning from start
                for i in 1..<len(global_data.xsects) {
                    if x > real(global_data.xsects[i].vertices[0].x) {
                        idx_loft_end = i+1
                    }
                }
                loft_end = real(global_data.xsects[idx_loft_end].vertices[0].x)
                create_cross_section_loft(&curr_loft, &global_data.xsects[idx_loft_end-1], &global_data.xsects[idx_loft_end])
            case .bbox:
                for i in 1..<len(global_data.bbox.corners) {
                    if x > real(global_data.bbox.corners[i].x) {
                        idx_loft_end = i+1
                    }
                }
                loft_end = real(global_data.bbox.corners[idx_loft_end].x)
                curr_rail = create_bbox_rail(&global_data.bbox, idx_loft_end)
            }
        }
        switch cfg.grid_parameterisation {
        case .rtheta:
            create_cross_section(&curr_xsect, &curr_loft, x)
            compute_grid_2d(&dn_grid, &up_grid, &global_data.rtheta_grid, &curr_xsect)
        case .bbox:
            t := bezier_t_from_x(curr_rail.bezier, x)
            corner := bezier_eval(curr_rail.bezier, t)
            compute_grid_2d_from_bbox(&dn_grid, &up_grid, &global_data.bbox_grid, corner)
        }
        add_3d_slice_of_hexes(&up_grid, &dn_grid)

        // Replace up_grid with dn_grid for next step.
        copy_grid_2d(&up_grid, &dn_grid)
    }

    if cfg.grid_parameterisation == .rtheta {
        delete_cross_section_loft(&curr_loft)
        delete_cross_section(&curr_xsect)
    }
    return true
}

