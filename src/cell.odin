package pingvin

Cell_Id :: distinct int

Cell :: struct {
    id :        Cell_Id,
    // Geometry
    volume:     f64,
    centroid:   Vector3,
    // Finite-volume elements
    faces :     [InterfaceLabels]Interface_Id,
    outsigns:   [InterfaceLabels]i8,
    // Flow and solver update elements
    cqs :       [Conserved_Quantities]complex128,
    pqs :       [Primitive_Quantities]complex128,
    R :         [Conserved_Quantities]complex128,
    dU :        [Conserved_Quantities]complex128,
}
