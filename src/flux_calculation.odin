package pingvin

import cmplx "core:math/cmplx"
import "core:fmt"

Flux_calculator :: enum {rusanov, van_leer}

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

// Delegator function to pick between Rusanov or van Leer fluxes
flux_calc :: proc (L, R: [Primitive_Quantities]complex128) -> (flux: [Conserved_Quantities]complex128) {
    switch globals.cfg.flux_calculator {
    case .rusanov:
        return rusanov_flux(L, R)
    case .van_leer:
        return van_leer_flux(L, R)
    }
    // compiler demands a return, so just use rusanov
    return rusanov_flux(L, R)
}

// Rusanov flux calculator.
// See pp 301--302 in Toro
rusanov_flux :: proc (L, R: [Primitive_Quantities]complex128) -> (flux: [Conserved_Quantities]complex128) {
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
        f = 0.5*(F_L[i] + F_R[i]) - 0.5*S_plus*(U_R[i] - U_L[i])
    }
    return flux    
}

// Using Toro's notatation
van_leer_flux :: proc (L, R: [Primitive_Quantities]complex128) -> (flux: [Conserved_Quantities]complex128) {
    // Compute f+ from state L
    rho := L[.rho]
    u := L[.xvel]
    v := L[.yvel]
    w := L[.zvel]
    gamma := globals.gamma
    a := cmplx.sqrt(gamma*L[.p]/rho)
    M := u/a

    f_plus : [Conserved_Quantities]complex128

    if real(M) >= 1.0 {
        f_plus = flux_vector(L)
    }
    else if real(M) <= -1.0 {
        // Do nothing, zero initialised
    }
    else {
        tmp := ((gamma - 1.)/2.)*M + 1.
        f_plus[.mass] = (1./4.)*rho*a*(1. + M)*(1. + M)
        f_plus[.xmom] = f_plus[.mass]*(2.*a/gamma)*tmp
        f_plus[.ymom] = f_plus[.mass]*v
        f_plus[.zmom] = f_plus[.mass]*w
        f_plus[.energy] = f_plus[.mass]*(((2.*a*a)/(gamma*gamma - 1.))*tmp*tmp + (1./2.)*(v*v + w*w))
    }

    // Compute f- from state R
    rho = R[.rho]
    u = R[.xvel]
    v = R[.yvel]
    w = R[.zvel]
    a = cmplx.sqrt(gamma*R[.p]/rho)
    M = u/a

    f_minus : [Conserved_Quantities]complex128

    if real(M) >= 1.0 {
        // Do nothing, since already zero initialised
    }
    else if real(M) <= -1.0 {
        f_minus = flux_vector(R)
    }
    else {
        tmp := ((gamma - 1.)/2.)*M - 1.
        f_minus[.mass] = -(1./4.)*rho*a*(1. - M)*(1. - M)
        f_minus[.xmom] = f_minus[.mass]*(2.*a/gamma)*tmp
        f_minus[.ymom] = f_minus[.mass]*v
        f_minus[.zmom] = f_minus[.mass]*w
        f_minus[.energy] = f_minus[.mass]*(((2.*a*a)/(gamma*gamma - 1.))*tmp*tmp + (1./2.)*(v*v + w*w))
    }

    for &f, i in flux {
        f = f_plus[i] + f_minus[i]
    }
    return flux
}
