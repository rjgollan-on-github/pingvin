package pingvin

import "core:math"
import "core:fmt"
import "core:os"

magnitude_yz :: proc (v: Vector3) -> f64 {
    return math.sqrt(v.z*v.z + v.y*v.y)
}

argument_yz :: proc (v: Vector3) -> f64 {
    return math.atan2(v.y, v.z)
}

find_cross_section_points :: proc (xsect: ^Cross_Section, theta: f64) -> (pA, pB: Vector3, result: bool) {
    for i in 0..<len(xsect.vertices)-1 {
        p0 := xsect.vertices[i]
        p1 := xsect.vertices[i+1]
        theta0 := math.atan2(p0.y, p0.z)
        theta1 := math.atan2(p1.y, p1.z)
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

compute_grid_2d :: proc (g: ^Grid_2d, rtg: ^Grid_rtheta, xsect: ^Cross_Section) {
    x := xsect.vertices[0].x
    for i in 0..<len(rtg.r_bar) {
        theta := rtg.theta[i]
        pA, pB, ok := find_cross_section_points(xsect, theta)
        if !ok {
            fmt.printfln("Error in compute_grid_2d: theta= %v", theta)
            os.exit(1)
        }
        theta0 := math.atan2(pA.y, pA.z)
        theta1 := math.atan2(pB.y, pB.z)
        w := (theta - theta0)/(theta1 - theta0)
        r_b := (1 - w)*magnitude_yz(pA) + w*magnitude_yz(pB)
        r := rtg.r_bar[i]*r_b
        z := r*math.cos(theta)
        y := r*math.sin(theta)
        append(&(global_data.vertices), Vector3{x, y, z})
        append(&(g.vertices), VtxInt(len(global_data.vertices)-1))
    }
}
