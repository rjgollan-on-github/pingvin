package pingvin

import cmplx "core:math/cmplx"

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
    flux[.mass] = -1.0*rho*u
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
        cell.pqs[.ke] = ke
        cell.pqs[.xvel] = u
    }
}

apply_downstream_flux :: proc (faces: #soa[]Interface) {
    for &f in faces {
        rho := f.right[.rho]
        p := f.right[.p]
        E := f.right[.e] + f.right[.ke]
        u := f.right[.xvel]
        v := f.right[.yvel]
        w := f.right[.zvel]
        f.flux[.mass] = rho*u
        f.flux[.xmom] = rho*u*u + p
        f.flux[.ymom] = rho*u*v
        f.flux[.zmom] = rho*u*w
        f.flux[.energy] = rho*E*u + p*u
    }
}

apply_slip_wall_flux :: proc (faces : []Interface_Id) {
    for f_id in faces {
        face := &global_data.x_faces[f_id]
        nx := complex(face.normal.x, 0.0)
        ny := complex(face.normal.y, 0.0)
        nz := complex(face.normal.z, 0.0)
        p := face.left[.p]
        face.flux[.mass] = complex(0.0, 0.0)
        face.flux[.energy] = complex(0.0, 0.0) 
        if face.left_cells[0] >= 0 {
            face.flux[.xmom] = -1.0*p*nx
            face.flux[.ymom] = -1.0*p*ny
            face.flux[.zmom] = -1.0*p*nz
        }
        else {
            face.flux[.xmom] = p*nx
            face.flux[.ymom] = p*ny
            face.flux[.zmom] = p*nz
        }
    }
}

apply_symm_flux :: apply_slip_wall_flux
