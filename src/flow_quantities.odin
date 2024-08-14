package pingvin

Conserved_Quantities :: enum {
    mass,
    xmom,
    ymom,
    zmom,
    energy,
}

cq :: Conserved_Quantities

Primitive_Quantities :: enum {
    rho, // density, kg/m^3
    p, // pressure, Pa
    T, // temperature, K
    e, // internal energy, J/kg
    ke, // kinetic energy, J/kg
    xvel, // x-velocity component, m/s
    yvel, // y-velocity component, m/s
    zvel, // z-velocity component, m/s
}

cq_from_prim :: proc (P : [Primitive_Quantities]complex128) -> (U: [Conserved_Quantities]complex128) {
    rho := P[.rho]
    u := P[.xvel]
    v := P[.yvel]
    w := P[.zvel]
    p := P[.p]
    E := P[.e] + P[.ke]
    U[.mass] = rho
    U[.xmom] = rho*u
    U[.ymom] = rho*v
    U[.zmom] = rho*w
    U[.energy] = E
    return U
}

prim_from_cq :: proc (U: [Conserved_Quantities]complex128) -> (P: [Primitive_Quantities]complex128) {
    rho := U[.mass]
    P[.rho] = rho
    rho_inv := 1.0/rho
    P[.xvel] = U[.xmom]*rho_inv
    P[.yvel] = U[.ymom]*rho_inv
    P[.zvel] = U[.zmom]*rho_inv
    P[.ke] = 0.5*(P[.xvel]*P[.xvel] + P[.yvel]*P[.yvel] + P[.zvel]*P[.zvel])
    E := U[.energy]*rho_inv
    P[.e] = E - P[.ke]
    P[.T] = P[.e]*(globals.gamma - 1)/globals.R_gas
    P[.p] = P[.rho]*globals.R_gas*P[.T]
    return P
}
