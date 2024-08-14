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

conserved_quantities_vector :: proc (P : [Primitive_Quantities]complex128) -> (U: [5]complex128) {
    rho := P[.rho]
    u := P[.xvel]
    v := P[.yvel]
    w := P[.zvel]
    p := P[.p]
    E := P[.e] + P[.ke]
    U[cq.mass] = rho
    U[cq.xmom] = rho*u
    U[cq.ymom] = rho*v
    U[cq.zmom] = rho*w
    U[cq.energy] = E
    return U
}


