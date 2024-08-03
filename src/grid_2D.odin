package pingvin

import "core:os"
import "core:strings"
import "core:fmt"
import "core:strconv"

Boundary :: enum{wall, symm}

Boundary_2d :: struct {
    faces: [dynamic]Face2,
}

Quad :: distinct [4]int

Grid_2d :: struct {
    vertices: [dynamic]Vector2,
    quads: [dynamic]Quad,
    wall_boundary: Boundary_2d,
    symm_boundary: Boundary_2d,
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

read_su2_2d_file :: proc (filepath: string, grid: ^Grid_2d) -> (err: GridInitError) {
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
    grid.vertices = make([dynamic]Vector2, n_points)
    for line in lines[line_points_start:][:n_points] {
        tokens = strings.fields(line)
        if len(tokens) != 3 {
            msg := fmt.tprintf("Error reading points in SU2 file '%s'. Three items expected on line.", filepath)
            return GridReadingError{msg=msg}
        }
        x := strconv.atof(tokens[0])
        y := strconv.atof(tokens[1])
        i := strconv.atoi(tokens[2])
        grid.vertices[i] = {x, y}
    }

    // 2. Find elements
    line_elem_start, n_elem : int
    line_elem_start, n_elem, ok = find_marker(lines, "NELEM")
    if !ok {
        msg := fmt.tprintf("Error reading NELEM in SU2 file '%s'", filepath)
        return GridReadingError{msg=msg}
    }
    grid.quads = make([dynamic]Quad, n_elem)
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
    grid.wall_boundary.faces = make([dynamic]Face2, n_wall_elems)
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
    grid.symm_boundary.faces = make([dynamic]Face2, n_symm_elems)
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



