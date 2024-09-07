package pingvin

import "core:fmt"

Cell_Id :: distinct int

Cell :: struct {
    id :        Cell_Id,
    // Geometry
    volume:     complex128,
    centroid:   Vector3,
    // Finite-volume elements
    faces :     [InterfaceLabels]Interface_Id,
    outsigns:   [InterfaceLabels]complex128,
    // Flow and solver update elements
    cqs :       [Conserved_Quantities]complex128,
    cqs0:       [Conserved_Quantities]f64,
    pqs :       [Primitive_Quantities]complex128,
    // Use 'R' to store current residual evaluation.
    // This may be a perturbed state.
    R :         [Conserved_Quantities]complex128,
    // Use 'Ru' to store unperturbed state.
    // Typically, we'll copy into this array once
    // at the start of each linear solve.
    Ru :        [Conserved_Quantities]f64,
    dU :        [Conserved_Quantities]f64,
}

compute_residual :: proc (cell: ^Cell) {
    sigma := 1.0e-250
    
    vol_inv := 1.0/cell.volume
    f_im := global_data.m_faces[cell.faces[.i_minus]]
    f_ip := global_data.m_faces[cell.faces[.i_plus]]
    f_jm := global_data.x_faces[cell.faces[.j_minus]]
    f_jp := global_data.x_faces[cell.faces[.j_plus]]
    f_km := global_data.x_faces[cell.faces[.k_minus]]
    f_kp := global_data.x_faces[cell.faces[.k_plus]]
    for cq in Conserved_Quantities {
        cell.R[cq] = f_im.flux[cq]*f_im.area*cell.outsigns[.i_minus]
        cell.R[cq] += f_ip.flux[cq]*f_ip.area*cell.outsigns[.i_plus]
        cell.R[cq] += f_jm.flux[cq]*f_jm.area*cell.outsigns[.j_minus]
        cell.R[cq] += f_jp.flux[cq]*f_jp.area*cell.outsigns[.j_plus]
        cell.R[cq] += f_km.flux[cq]*f_km.area*cell.outsigns[.k_minus]
        cell.R[cq] += f_kp.flux[cq]*f_kp.area*cell.outsigns[.k_plus]
        
        cell.R[cq] *= -1.0*vol_inv
    }
}

copy_residual_to_base_state :: proc (cell: ^Cell) {
    for cq in Conserved_Quantities do cell.Ru[cq] = real(cell.R[cq])
    return
}

copy_cqs_to_base_state :: proc (cell: ^Cell) {
    for cq in Conserved_Quantities do cell.cqs0[cq] = real(cell.cqs[cq])
    return
}
