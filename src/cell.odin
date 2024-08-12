package pingvin

Cell_Id :: distinct int

Cell :: struct {
    id :        Cell_Id,
    cqs :       [Conserved_Quantities]complex128,
    pqs :       [Primitive_Quantities]complex128,
    faces :     [InterfaceLabels]Interface_Id,
    outsigns:   [InterfaceLabels]i8,
    R :         [Conserved_Quantities]complex128,
    dU :        [Conserved_Quantities]complex128,
}
