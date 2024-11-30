package pingvin

import "core:fmt"

Interface_Id :: distinct int

InterfaceLabels :: enum {i_minus, i_plus, j_minus, j_plus, k_minus, k_plus}

InterfaceType :: enum {interior, wall, symm, off_wall}

Interface :: struct {
    area :        complex128,
    normal :      Vector3,
    t1 :          Vector3,
    t2 :          Vector3,
    centroid :    Vector3,
    flux :        [Conserved_Quantities]complex128,
    left :        [Primitive_Quantities]complex128,
    right :       [Primitive_Quantities]complex128,
    left_cells :  [2]Cell_Id,
    right_cells : [2]Cell_Id,
    quad :        Quad,
}

transform_interior_states_to_local_frame :: proc (face_ids: []Interface_Id) {
    for f_id in face_ids {
        f := &global_data.x_faces[f_id]
        transform_pq_vel_to_local_frame(&f.left, f.normal, f.t1, f.t2)
        transform_pq_vel_to_local_frame(&f.right, f.normal, f.t1, f.t2)
    }
}

transform_interior_flux_to_global_frame :: proc (face_ids: []Interface_Id) {
    for f_id in face_ids {
        f := &global_data.x_faces[f_id]
        transform_flux_to_global_frame(&f.flux, f.normal, f.t1, f.t2)
    }
}
