package pingvin

import "core:fmt"
import "core:os"
import "core:strings"
import "core:strconv"
import "core:log"


Cross_Section :: struct {
    vertices: [dynamic]Vector3,
    slopes: [dynamic]Vector3,
}

CrossSectionError :: union {
    FileOpenFailed, 
}

read_cross_section :: proc (xsect: ^Cross_Section, dir: string, isect: int) -> (err: CrossSectionError) {
    filename := fmt.tprintf("%s/xsect-%03d", dir, isect)
    data, ok := os.read_entire_file(filename)
    log.debugf("filename= %v ok= %v", filename, ok)
    if !ok {
        return FileOpenFailed{filename=filename}
    }
    defer delete(data)

    lines := strings.split_lines(string(data))

    // 0. First line should say how many entries
    tokens := strings.fields(lines[0])
    n_entries := strconv.atoi(tokens[0])

    // 1. Now gather vertices and slopes
    xsect.vertices = make([dynamic]Vector3, n_entries)
    xsect.slopes = make([dynamic]Vector3, n_entries)
    // Line format expected is:
    // x  y  z  dydx  dydz
    for i in 1..=n_entries {
        tokens = strings.fields(lines[i])
        xsect.vertices[i-1] = {strconv.atof(tokens[0]),
                               strconv.atof(tokens[1]),
                               strconv.atof(tokens[2])}
        xsect.slopes[i-1] = {0.0,
                             strconv.atof(tokens[3]),
                             strconv.atof(tokens[4])}
    }

    return nil
}


