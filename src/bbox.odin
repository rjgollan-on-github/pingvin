package pingvin

import "core:fmt"
import "core:os"
import "core:strings"
import "core:strconv"

Bounding_box :: struct {
    corners : []Vector3,
    slopes : []Vector3,
}

Bbox_rail :: struct {
    bezier : CubicBezier,
    start : Vector3,
    end : Vector3,
}

read_bbox :: proc (filename : string) {
    data, ok := os.read_entire_file(filename)
    if !ok {
        fmt.printfln("pvn: Failed to read file: %s", filename)
        fmt.printfln("pvn: Exiting.")
        os.exit(1)
    }
    defer delete(data)

    lines := strings.split_lines(string(data))

    // 0. First line should say how many entries
    tokens := strings.fields(lines[0])
    n_entries := strconv.atoi(tokens[0])

    global_data.bbox.corners = make([]Vector3, n_entries)
    global_data.bbox.slopes = make([]Vector3, n_entries)

    // 1. Read corners and slopes
    // Line format expected is:
    // x  y  z  dydx  dydz
    for i in 1..=n_entries {
        tokens = strings.fields(lines[i])
        global_data.bbox.corners[i-1] = {strconv.atof(tokens[0]),
                                         strconv.atof(tokens[1]),
                                         strconv.atof(tokens[2])}
        global_data.bbox.slopes[i-1] = {1.0,
                                        strconv.atof(tokens[3]),
                                        strconv.atof(tokens[4])}
    }

    return
}

create_bbox_rail :: proc (rail : ^Bbox_rail, bbox : ^Bounding_box, upper_idx : int) {
    i := upper_idx - 1
    b0 := bbox.corners[i]
    b3 := bbox.corners[i+1]
    dx := b3.x - b0.x
    b1 := b0 + (1./3)*dx*bbox.slopes[i]
    b2 := b3 - (1./3)*dx*bbox.slopes[i+1]
    rail.bezier = CubicBezier{b0, b1, b2, b3}
    rail.start = b0
    rail.end = b3
}

update_rail :: proc (rail: ^Bbox_rail, x: f64) {
    idx_rail_end : int
    if x > real(rail.end.x) {
        for i in 1..<len(global_data.bbox.corners)-1 {
            if x > real(global_data.bbox.corners[i].x) {
                idx_rail_end = i + 1
            }
        }
        create_bbox_rail(rail, &global_data.bbox, idx_rail_end)
    }
    return
}
