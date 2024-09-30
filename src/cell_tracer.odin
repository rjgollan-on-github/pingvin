package pingvin

import "core:fmt"
import "core:os"

Cell_tracer :: struct {
    vols :  [dynamic]f64,
    ctrs :  [dynamic]Vector3,
    pqs  :  [dynamic][Primitive_Quantities]f64,
    Ms   :  [dynamic]f64,
}


write_cell_trace :: proc (filename : string, ct : Cell_tracer) {
    f, err := os.open(filename, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o644)
    defer os.close(f)

    fmt.fprintln(f, "x  y  z  vol  rho  p  T  e  xvel  yvel  zvel  Mach")

    for i in 0..<len(ct.vols) {
        fmt.fprintf(f, "%.6e %.6e %.6e %.6e ", real(ct.ctrs[i].x), real(ct.ctrs[i].y), real(ct.ctrs[i].z), ct.vols[i])
        for q in Primitive_Quantities {
            fmt.fprintf(f, "%.18e ", ct.pqs[i][q])
        }
        fmt.fprintfln(f, "%.18e", ct.Ms[i])
    }
    return    
}


