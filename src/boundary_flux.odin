package pingvin

import cmplx "core:math/cmplx"
import "core:fmt"

apply_inflow :: proc (slice: ^Slice) {
    cfg := globals.cfg
    p := complex(cfg.p_inflow, 0)
    T := complex(cfg.T_inflow, 0)
    M := complex(cfg.Mach_inflow, 0)
    gamma := globals.gamma
    R_gas := globals.R_gas

    rho := p/(R_gas*T)
    a := cmplx.sqrt(gamma*R_gas*T)
    u := M*a
    e := (R_gas/(gamma - 1.0))*T
    ke := 0.5*u*u
    E := e + ke
    flux : [Conserved_Quantities]complex128
    flux[.mass] = -rho*u
    flux[.xmom] = flux[.mass]*u - p
    flux[.energy] = flux[.mass]*E - p*u
    for &f in slice.up_faces {
        f.flux = flux
    }
    for &cell in global_data.cells[slice.first_cell:slice.last_cell+1] {
        cell.pqs[.rho] = rho
        cell.pqs[.p] = p
        cell.pqs[.T] = T
        cell.pqs[.e] = e
        cell.pqs[.xvel] = u
        cell.cqs = cq_from_prim(cell.pqs)
    }
}

apply_downstream_flux :: proc (faces: #soa[]Interface) {
    for &f in faces {
        rho := f.right[.rho]
        p := f.right[.p]
        u := f.right[.xvel]
        v := f.right[.yvel]
        w := f.right[.zvel]
        ke := 0.5*(u*u + v*v + w*w)
        E := f.right[.e] + ke
        f.flux[.mass] = -rho*u
        f.flux[.xmom] = f.flux[.mass]*u - p
        f.flux[.ymom] = f.flux[.mass]*v
        f.flux[.zmom] = f.flux[.mass]*w
        f.flux[.energy] = f.flux[.mass]*E - p*u
    }
}

bap_downstream_flux :: proc (P, Pm1 : [Primitive_Quantities]complex128) -> (flux : [Conserved_Quantities]complex128) {
    dx := complex(globals.cfg.dx, 0.0)
    // mass
    f_j_mass := -P[.rho] * P[.xvel]
    f_jm1_mass := -Pm1[.rho] * Pm1[.xvel]
    s_l := (f_j_mass - f_jm1_mass)/dx
    s_r := complex(0.0, 0.0)
    s_j := slope_average(s_l, s_r)
    flux[.mass] = f_j_mass + s_l*dx/2.0
    // x-mom
    f_j := f_j_mass * P[.xvel] - P[.p]
    f_jm1 := f_jm1_mass * Pm1[.xvel] - Pm1[.p]
    s_l = (f_j - f_jm1)/dx
    s_j = slope_average(s_l, s_r)
    flux[.xmom] = f_j + s_l*dx/2.0
    // y-mom
    f_j = f_j_mass * P[.yvel]
    f_jm1 = f_jm1_mass * Pm1[.yvel]
    s_l = (f_j - f_jm1)/dx
    s_j = slope_average(s_l, s_r)
    flux[.ymom] = f_j + s_l*dx/2.0
    // z-mom
    f_j = f_j_mass * P[.zvel]
    f_jm1 = f_jm1_mass * Pm1[.zvel]
    s_l = (f_j - f_jm1)/dx
    s_j = slope_average(s_l, s_r)
    flux[.zmom] = f_j + s_l*dx/2.0
    // energy
    E := P[.e] + 0.5*(P[.xvel]*P[.xvel] + P[.yvel]*P[.yvel] + P[.zvel]*P[.zvel])
    f_j = f_j_mass * E - P[.p]*P[.xvel]
    E = Pm1[.e] + 0.5*(Pm1[.xvel]*Pm1[.xvel] + Pm1[.yvel]*Pm1[.yvel] + Pm1[.zvel]*Pm1[.zvel])
    f_jm1 = f_j_mass * E - Pm1[.p]*Pm1[.xvel]
    s_l = (f_j - f_jm1)/dx
    s_j = slope_average(s_l, s_r)
    flux[.energy] = f_j + s_l*dx/2.0

    return flux
}

apply_slip_wall_flux :: proc (faces : []Interface_Id) {
    for f_id in faces {
        face := &global_data.x_faces[f_id]
        nx := face.normal.x
        ny := face.normal.y
        nz := face.normal.z
        face.flux[.mass] = complex(0.0, 0.0)
        face.flux[.energy] = complex(0.0, 0.0) 
        if face.left_cells[0] >= 0 {
            p := face.left[.p]
            face.flux[.xmom] = p*nx
            face.flux[.ymom] = p*ny
            face.flux[.zmom] = p*nz
        }
        else {
            p := face.right[.p]
            face.flux[.xmom] = p*nx
            face.flux[.ymom] = p*ny
            face.flux[.zmom] = p*nz
        }
    }
}

apply_symm_flux :: apply_slip_wall_flux
