package pingvin

import cmplx "core:math/cmplx"
import "core:fmt"

abs_c128 :: proc(z : complex128) -> complex128 {
    if real(z) >= 0.0 {
        return z
    }
    else {
        return complex(-real(z), -imag(z))
    }
}


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
    ke := 0.5*(P[.xvel]*P[.xvel] + P[.yvel]*P[.yvel] + P[.zvel]*P[.zvel])
    p := P[.p]
    E := P[.e] + ke
    flux[.mass] = rho*u
    flux[.xmom] = rho*u*u + p
    flux[.ymom] = rho*u*v
    flux[.zmom] = rho*u*w
    flux[.energy] = rho*u*E + p*u
    return flux
}

// Rusanov flux calculator.
// See pp 301--302 in Toro
flux_calc :: proc (L, R: [Primitive_Quantities]complex128) -> (flux: [Conserved_Quantities]complex128) {
    F_L := flux_vector(L)
    F_R := flux_vector(R)
    U_L := cq_from_prim(L)
    U_R := cq_from_prim(R)

    u_L := L[.xvel]
    u_R := R[.xvel]
    gamma := globals.gamma
    a_L := cmplx.sqrt(gamma*L[.p]/L[.rho]) 
    a_R := cmplx.sqrt(gamma*R[.p]/R[.rho])
    S_plus := max_complex(abs_c128(u_L - a_L), abs_c128(u_R - a_R), abs_c128(u_L + a_L), abs_c128(u_R + a_R))

    for &f, i in flux {
        f = 0.5*(F_L[i] + F_R[i]) //- 0.5*S_plus*(U_R[i] - U_L[i])
    }
    return flux    
}
