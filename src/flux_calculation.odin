package pingvin

import cmplx "core:math/cmplx"

max_c128 :: proc(z1, z2: complex128) -> complex128 {
    if real(z1) >= real(z2) {
        return z1
    }
    else {
        return z2
    }
}
max_complex :: proc(values: ..complex128) -> complex128 {
    mv := values[0]
    for v in values[1:] {
        mv = max_c128(mv, v)
    }
    return mv
}

flux_vector :: proc(P: [Primitive_Quantities]complex128) -> (flux: [Conserved_Quantities]complex128) {
    rho := P[.rho]
    u := P[.xvel]
    v := P[.yvel]
    w := P[.zvel]
    p := P[.p]
    E := P[.e] + P[.ke]
    flux[.mass] = rho*u
    flux[.xmom] = rho*u*u + p
    flux[.ymom] = rho*u*v
    flux[.zmom] = rho*u*w
    flux[.energy] = u*(E + p)
    return flux
}

// Rusanov flux calculator.
// See pp 301--302 in Toro
flux_calc :: proc (L, R: [Primitive_Quantities]complex128) -> (flux: [Conserved_Quantities]complex128) {
    F_L := flux_vector(L)
    F_R := flux_vector(R)
    U_L := conserved_quantities_vector(L)
    U_R := conserved_quantities_vector(R)

    u_L := L[.xvel]
    u_R := R[.xvel]
    gamma := globals.gamma
    R_gas := globals.R_gas
    a_L := cmplx.sqrt(gamma*R_gas*L[.p]/L[.rho]) 
    a_R := cmplx.sqrt(gamma*R_gas*R[.p]/R[.rho])

    S_plus := max_complex(abs(u_L - a_L), abs(u_R - a_R), abs(u_L + a_L), abs(u_R + a_R))

    for &f, i in flux {
        f = 0.5*(F_L[i] + F_R[i]) - 0.5*S_plus*(U_L[i] - U_R[i])
    }
    return flux    
}
