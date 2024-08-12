package pingvin

import "core:testing"
import "core:log"


tet_dipyramid_volume :: proc(p0, p1, p2, p3, pc: Vector3) -> (volume: f64) {
    pb := 0.25*(p0 + p1 + p2 + p3)
    A := pc - pb
    B := p0 - p1 + p3 - p2
    C := p0 + p1 - p3 - p2
    volume = A.x*(B.y*C.z - B.z*C.y)
    volume += -B.x*(A.y*C.z - A.z*C.y)
    volume += C.x*(A.y*B.z - A.z*C.y)
    return volume
}

@(test)
cell_formation_test :: proc (t: ^testing.T) {
    global_data.vertices = make([dynamic]Vector3, 8)
    defer delete(global_data.vertices)
    v := global_data.vertices[:]
    v[0] = {0.0000000000000000e+00, 2.9826897612595271e-01, 4.6100533824970175e-01}
    v[1] = {0.0000000000000000e+00, 2.9364899967529617e-01, 6.1482333830801961e-01}
    v[2] = {0.0000000000000000e+00, 4.3312491699843125e-01, 5.8658766497157921e-01}
    v[3] = {0.0000000000000000e+00, 4.4386251683769462e-01, 4.4386254493036931e-01}
    v[4] = {6.4399999999999979e-01, 2.9743075457556650e-01, 4.5970978074860647e-01}
    v[5] = {6.4399999999999979e-01, 2.9282376159801576e-01, 6.1309550804081803e-01}
    v[6] = {6.4399999999999979e-01, 4.3190771098736402e-01, 5.8493918435132464e-01}
    v[7] = {6.4399999999999979e-01, 4.4261512902442718e-01, 4.4261515703815296e-01}

    pc := (1./8)*(v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7])

    // Check orientation for i_minus
    vol := tet_dipyramid_volume(v[0], v[1], v[2], v[3], pc)
    testing.expect(t, vol < 0.0)

    // i_plus
    vol = tet_dipyramid_volume(v[4], v[5], v[6], v[7], pc)
    testing.expect(t, vol > 0.0)

    // j_minus
    vol = tet_dipyramid_volume(v[0], v[4], v[5], v[1], pc)
    testing.expect(t, vol < 0.0)

    // j_plus
    vol = tet_dipyramid_volume(v[2], v[6], v[7], v[3], pc)
    testing.expect(t, vol < 0.0)

    // k_minus
    vol = tet_dipyramid_volume(v[3], v[7], v[4], v[0], pc)
    testing.expect(t, vol < 0.0)

    // k_plus
    vol = tet_dipyramid_volume(v[1], v[5], v[6], v[2], pc)
    testing.expect(t, vol < 0.0)
    
}    
