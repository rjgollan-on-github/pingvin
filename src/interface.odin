package pingvin

Interface_Id :: distinct int

InterfaceLabels :: enum {i_plus, i_minus, j_plus, j_minus, k_plus, k_minus}

Interface :: struct {
    area :        f64,
    normal :      Vector3,
    flux :        [Conserved_Quantities]complex128,
    left :        [Primitive_Quantities]complex128,
    right :       [Primitive_Quantities]complex128,
    left_cells :  [2]Cell_Id,
    right_cells : [2]Cell_Id,
}
