package pingvin

Cell_Id :: distinct int

Cell :: struct {
    id :        Cell_Id,
    // Geometry
    volume:     f64,
    centroid:   Vector3,
    // Finite-volume elements
    faces :     [InterfaceLabels]Interface_Id,
    outsigns:   [InterfaceLabels]complex128,
    // Flow and solver update elements
    cqs :       [Conserved_Quantities]complex128,
    pqs :       [Primitive_Quantities]complex128,
    R :         [Conserved_Quantities]complex128,
    dU :        [Conserved_Quantities]complex128,
}

compute_residual :: proc (cell: ^Cell) {
    vol_inv := 1.0/complex(cell.volume, 0.0)
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
