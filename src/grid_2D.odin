package pingvin

import "core:os"
import "core:strings"
import "core:fmt"
import "core:strconv"
import "core:flags"
import "core:log"
import "core:math"

Boundary :: enum{wall, symm}

Boundary_2d :: struct {
    faces: [dynamic]Face2,
}

Quad :: distinct [4]int

Grid_2d :: struct {
    vertices: [dynamic]VtxInt,
    quads: []Quad,
    wall_boundary: ^Boundary_2d,
    symm_boundary: ^Boundary_2d,
}

Grid_rtheta :: struct {
    r_bar: [dynamic]f64,
    theta: [dynamic]f64,
}

GridInitError :: union {
    FileOpenFailed,
    GridReadingError,
    GridFormationError,
}

FileOpenFailed :: struct {
    filename: string,
}

GridReadingError :: struct {
    msg: string,
}

GridFormationError :: struct {
    msg: string,
}

GenGridCmd := Command {
    main = generate_grid,
    description = "Generate 3D grid from 2D grid and cross sections.",
    short_desc = "...",
}

GenGridOptions :: struct {
    job_file: string `args:"name=job" usage:"job script (Lua file) default: job.lua"`,
}

generate_grid :: proc (args: []string) -> (result: bool) {
    // 0. Some command line and config setup
    cmd_name := args[0]
    opt := GenGridOptions{job_file="job.lua"}
    err := flags.parse(&opt, args[1:])

    if err != nil {
        flags.print_errors(typeid_of(GenGridOptions), err, cmd_name)
        os.exit(1)
    }
    grid : Grid_2d
    cfg := read_config_from_lua_file(opt.job_file)

    // 1. Read in grid at initial plane
    read_su2_2d_file(&grid, cfg.grid2d_file)
    log.debugf("Number of points= %d", len(grid.vertices))
    log.debugf("Number of quads= %d", len(grid.quads))
    log.debugf("Number of wall elems= %d", len(grid.wall_boundary.faces))
    log.debugf("Numer of symm elems= %d", len(grid.symm_boundary.faces))

    // 2. Read in cross section at initial plane
    xsect : Cross_Section
    xsect_err := read_cross_section(&xsect, cfg.cross_section_dir, 0)
    switch subtype in xsect_err {
    case FileOpenFailed:
        fmt.printfln("Unable to open file with initial-plane cross section.")
        os.exit(1)
    }

    // 2. Prepare (global) r-theta grid
    allocate_rtheta_grid(len(grid.vertices))
    compute_rtheta_grid(&(global_data.rtheta_grid), &grid, &xsect)

    // 3. Read in next cross-section
    xsect1 : Cross_Section
    xsect_err = read_cross_section(&xsect1, cfg.cross_section_dir, 1)
    grid1 : Grid_2d
    grid1.quads = global_data.quads[:]
    grid1.wall_boundary = &(global_data.wall_boundary)
    grid1.symm_boundary = &(global_data.wall_boundary)
    compute_grid_2d(&grid1, &(global_data.rtheta_grid), &(xsect1))
    write_grid_2d_as_vtk("test-g1.vtk", &grid1)
    
    
    return true
}

find_marker :: proc (lines: []string, marker: string) -> (line_start, n_things: int, result: bool) {
    for line, line_no in lines {
        if strings.has_prefix(line, marker) {
            line_start = line_no + 1
            tokens := strings.fields(line)
            n_things = strconv.atoi(tokens[1])
            return line_start, n_things, true
        }
    }
    // If we get here, error
    return line_start, n_things, false
}

find_marker_tag :: proc (lines: []string, tag: string) -> (line_start, n_elems: int, result: bool) {
    match := fmt.tprintf("MARKER_TAG= %s", tag)
    for line, line_no in lines {
        if strings.compare(line, match) == 0 {
            line_start = line_no + 2
            tokens := strings.fields(lines[line_no+1])
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

    lines := strings.split_lines(string(data))
    // 0. Examine dimensions on first line
    tokens := strings.fields(lines[0])
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
        tokens = strings.fields(line)
        if len(tokens) != 3 {
            msg := fmt.tprintf("Error reading points in SU2 file '%s'. Three items expected on line.", filepath)
            return GridReadingError{msg=msg}
        }
        z := strconv.atof(tokens[0])
        y := strconv.atof(tokens[1])
        i := strconv.atoi(tokens[2])
        global_data.vertices[i] = {0.0, y, z}
        append(&(grid.vertices), VtxInt(i))
    }

    // 2. Find elements
    line_elem_start, n_elem : int
    line_elem_start, n_elem, ok = find_marker(lines, "NELEM")
    if !ok {
        msg := fmt.tprintf("Error reading NELEM in SU2 file '%s'", filepath)
        return GridReadingError{msg=msg}
    }
    allocate_quads(n_elem)
    grid.quads = global_data.quads[:]
    for line in lines[line_elem_start:][:n_elem] {
        tokens = strings.fields(line)
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
        grid.quads[i] = {strconv.atoi(tokens[1]),
                         strconv.atoi(tokens[2]),
                         strconv.atoi(tokens[3]),
                         strconv.atoi(tokens[4])}
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
        tokens = strings.fields(line)
        if len(tokens) != 3 {
            msg := fmt.tprintf("Error reading quad in SU2 file '%s'. Three items expected on line.", filepath)
            return GridReadingError{msg=msg}
        }
        e_type := strconv.atoi(tokens[0])
        if (e_type != int(VTKElement.line)) {
            msg := fmt.tprintf("Error in SU2 file '%s'. Line element expected (3); found %d.", filepath, e_type)
            return GridReadingError{msg=msg}
        }
        grid.wall_boundary.faces[i] = {strconv.atoi(tokens[1]), strconv.atoi(tokens[2])}
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
        tokens = strings.fields(line)
        if len(tokens) != 3 {
            msg := fmt.tprintf("Error reading quad in SU2 file '%s'. Three items expected on line.", filepath)
            return GridReadingError{msg=msg}
        }
        e_type := strconv.atoi(tokens[0])
        if (e_type != int(VTKElement.line)) {
            msg := fmt.tprintf("Error in SU2 file '%s'. Line element expected (3); found %d.", filepath, e_type)
            return GridReadingError{msg=msg}
        }
        grid.symm_boundary.faces[i] = {strconv.atoi(tokens[1]), strconv.atoi(tokens[2])}
    }

    return nil
}

uoflowz :: proc (q: f64, tiny : f64 = math.F32_MIN, huge : f64 = math.F32_MAX) -> (qsafe: f64) {
    qsafe = q
    if math.abs(q) < tiny { qsafe = 0.0 }
    if math.abs(q) > huge { qsafe = math.copy_sign(huge, q)}
    return q
}

write_grid_2d_as_vtk :: proc (filename: string, g: ^Grid_2d) {
    f, err := os.open(filename, os.O_WRONLY | os.O_CREATE, 0o644)
    defer os.close(f)

    fmt.fprintln(f, "# vtk DataFile Version 2.0")
    fmt.fprintln(f, "grid_2d, created by pingvin 0.0.1")
    fmt.fprintln(f, "ASCII")
    fmt.fprintln(f, "DATASET UNSTRUCTURED_GRID")
    fmt.fprintfln(f, "POINTS %d double", len(g.vertices))
    for v_id in g.vertices {
        v := global_data.vertices[v_id]
        fmt.fprintfln(f, "%.16e %.16e %.16e", uoflowz(v.x), uoflowz(v.y), uoflowz(v.z))
    }
    n_cells := len(g.quads)
    fmt.fprintln(f, "")
    fmt.fprintfln(f, "CELLS %d %d", n_cells, n_cells*5)
    for q in g.quads {
        fmt.fprintfln(f, "%d %d %d %d %d", 4, q[0], q[1], q[2], q[3])
    }
    fmt.fprintfln(f, "CELL_TYPES %d", n_cells)
    for i in 0..<n_cells {
        fmt.fprintfln(f, "%d", VTKElement.quad)
    }
    
}


