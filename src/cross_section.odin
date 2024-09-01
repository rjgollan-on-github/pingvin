package pingvin

import "core:fmt"
import "core:os"
import "core:strings"
import "core:strconv"
import "core:log"

TYPICAL_NO_CROSS_SECTIONS := 50

Cross_Section :: struct {
    vertices: [dynamic]Vector3,
    slopes: [dynamic]Vector3,
}

Cross_Section_Loft :: struct {
    end:     f64,
    beziers: [dynamic]CubicBezier,
}

CrossSectionError :: union {
    FileOpenFailed,
    DirectoryOpenFailed,
}

read_cross_section :: proc (xsect: ^Cross_Section, filename: string, isect: int) -> (err: CrossSectionError) {
    data, ok := os.read_entire_file(filename)
    if !ok {
        return FileOpenFailed{filename=filename}
    }
    defer delete(data)

    lines := strings.split_lines(string(data))

    // 0. First line should say how many entries
    tokens := strings.fields(lines[0])
    n_entries := strconv.atoi(tokens[0])

    // 1. Now gather vertices and slopes
    allocate_cross_section(xsect, n_entries)
    // Line format expected is:
    // x  y  z  dydx  dydz
    for i in 1..=n_entries {
        tokens = strings.fields(lines[i])
        xsect.vertices[i-1] = {strconv.atof(tokens[0]),
                               strconv.atof(tokens[1]),
                               strconv.atof(tokens[2])}
        xsect.slopes[i-1] = {1.0,
                             strconv.atof(tokens[3]),
                             strconv.atof(tokens[4])}
    }

    return nil
}

read_all_cross_sections :: proc (dir: string, n_xsects: int) -> (err: CrossSectionError) {
    global_data.xsects = make([dynamic]Cross_Section, n_xsects)
    for i in 0..<n_xsects {
        filename := fmt.tprintf("%s/xsect-%03d", dir, i)
        read_cross_section(&global_data.xsects[i], filename, i)
    }
    return nil
}

allocate_cross_section :: proc (xsect: ^Cross_Section, n: int) {
    xsect.vertices = make([dynamic]Vector3, n)
    xsect.slopes = make([dynamic]Vector3, n)
}

delete_cross_section :: proc(xsect: ^Cross_Section) {
    delete(xsect.vertices)
    delete(xsect.slopes)
}

allocate_cross_section_loft :: proc (loft: ^Cross_Section_Loft, n: int) {
    loft.beziers = make([dynamic]CubicBezier, n)
}

delete_cross_section_loft :: proc(loft: ^Cross_Section_Loft) {
    delete(loft.beziers)
}

update_loft :: proc(loft: ^Cross_Section_Loft, x: f64) {
    idx_loft_end : int
    if x > loft.end {
        // Search for new loft end in cross sections
        for i in 1..<len(global_data.xsects)-1 {
            if x > global_data.xsects[i].vertices[0].x {
                idx_loft_end = i + 1
            }
        }
        create_cross_section_loft(loft, &global_data.xsects[idx_loft_end-1], &global_data.xsects[idx_loft_end])
    }
}

create_cross_section_loft :: proc (loft: ^Cross_Section_Loft, up, dn: ^Cross_Section) {
    n_seg := len(up.vertices)
    for i in 0..<n_seg {
        b0 := up.vertices[i]
        b3 := dn.vertices[i]
        dx := b3.x - b0.x
        b1 := b0 + (1./3)*dx*up.slopes[i]
        b2 := b3 - (1./3)*dx*dn.slopes[i]
        loft.beziers[i] = CubicBezier{b0, b1, b2, b3}
    }
    loft.end = dn.vertices[0].x
}



create_cross_section :: proc (xsect: ^Cross_Section, loft: ^Cross_Section_Loft, x: f64) {
    for i in 0..<len(xsect.vertices) {
        t := bezier_t_from_x(loft.beziers[i], x)
        xsect.vertices[i] = bezier_eval(loft.beziers[i], t)
    }
}



