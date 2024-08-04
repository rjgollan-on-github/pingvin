package pingvin

import "core:math"
import "core:fmt"
import "core:os"

magnitude :: proc (v: Vector2) -> f64 {
    return math.sqrt(v.x*v.x + v.y*v.y)
}

magnitude_yz :: proc (v: Vector3) -> f64 {
    return math.sqrt(v.z*v.z + v.y*v.y)
}

argument :: proc (v: Vector2) -> f64 {
    return math.atan2(v.y, v.x)
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
    for v, i in g.vertices {
        r := magnitude(v)
        theta := argument(v)
        pA, pB, ok := find_cross_section_points(xsect, theta)
        if !ok {
            fmt.printfln("Error in compute_rtheta_grid: r= %v theta= %v", r, theta)
            os.exit(1)
        }
        theta0 := math.atan2(pA.y, pA.z)
        theta1 := math.atan2(pB.y, pB.z)
        w := (theta - theta0)/(theta1 - theta0)
        r_b := (1 - w)*magnitude_yz(pA) + w*magnitude_yz(pB)
        // Now populate the rtheta grid
        rtg.r_bar[i] = r/r_b
        rtg.theta[i] = theta
    }
}

compute_grid_2d :: proc (g: ^Grid_2d, rtg: ^Grid_rtheta, xsect: ^Cross_Section) {
    g.vertices = make([dynamic]Vector2, len(rtg.r_bar))
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
        x := r*math.cos(theta)
        y := r*math.sin(theta)
        g.vertices[i] = {x, y}
    }
}
