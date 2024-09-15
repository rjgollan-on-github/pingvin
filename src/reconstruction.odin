package pingvin

import "core:fmt"
import "core:math/cmplx"

// For interior faces
low_order_reconstruction :: proc (faces: []Interface_Id) {
    for f_id in faces {
        face := &global_data.x_faces[f_id]
        face.left = global_data.cells[face.left_cells[0]].pqs
        face.right = global_data.cells[face.right_cells[0]].pqs
    }
}

// For downstream faces
low_order_recon_downstream :: proc (faces: #soa[]Interface) {
    for &face in faces {
        face.right = global_data.cells[face.right_cells[0]].pqs
    }
}

// At slice bouondaries
low_order_recon_boundary :: proc (faces: []Interface_Id) {
    for f_id in faces {
        face := &global_data.x_faces[f_id]
        if face.left_cells[0] >= 0 {
            face.left = global_data.cells[face.left_cells[0]].pqs
        }
        else { // must be domain cell on right
            face.right = global_data.cells[face.right_cells[0]].pqs
        }
    }
}

// flux reconstruction via biased averaging procedure
B :: proc (x: complex128) -> complex128 {
    return x/(cmplx.sqrt(1. + x*x) + 1.)
}

B_inv :: proc (x: complex128) -> complex128 {
    return (2.*x)/(1. - x*x)
}

slope_average :: proc (s_l, s_r : complex128) -> complex128 {
    return B_inv((B(s_l) + B(s_r))/2.)
}
