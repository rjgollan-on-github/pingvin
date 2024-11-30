package pingvin

import "core:fmt"
import "core:math/cmplx"

// For interior faces
low_order_reconstruction :: proc (faces: []Interface_Id) {
    for f_id in faces {
        face := &global_data.x_faces[f_id]
        face.left = global_data.cells[face.left_cells[0]].pqs
        face.right = global_data.cells[face.right_cells[0]].pqs
    }
}

// For downstream faces
low_order_recon_downstream :: proc (faces: #soa[]Interface) {
    for &face in faces {
        face.right = global_data.cells[face.right_cells[0]].pqs
    }
}

// At slice bouondaries
low_order_recon_boundary :: proc (faces: []Interface_Id) {
    for f_id in faces {
        face := &global_data.x_faces[f_id]
        if face.left_cells[0] >= 0 {
            face.left = global_data.cells[face.left_cells[0]].pqs
        }
        else { // must be domain cell on right
            face.right = global_data.cells[face.right_cells[0]].pqs
        }
    }
}

high_order_recon_interior :: proc (faces: []Interface_Id) {
    for f_id in faces {
        face := &global_data.x_faces[f_id]
        L1 := face.left_cells[1]
        L0 := face.left_cells[0]
        R0 := face.right_cells[0]
        R1 := face.right_cells[1]
        //fmt.printfln("f_id= %d: L1= %d L0= %d R0= %d R1= %d", f_id, L1, L0, R0, R1)
        dx0 := magnitude(global_data.cells[L0].centroid - global_data.cells[L1].centroid)
        dx1 := magnitude(global_data.cells[R0].centroid - global_data.cells[L0].centroid)
        dx2 := magnitude(global_data.cells[R1].centroid - global_data.cells[R0].centroid)
        dxL := magnitude(face.centroid - global_data.cells[L0].centroid)
        dxR := magnitude(face.centroid - global_data.cells[R0].centroid)
        //fmt.printfln("dx0= %v dx1= %v dx2= %v dxL= %v dxR= %v", dx0, dx1, dx2, dxL, dxR)
        // Reconstruct rho, p and velocities
        s0, s1, s2 : complex128
        sL, sR : complex128
        // rho
        s0 = (global_data.cells[L0].pqs[.rho] - global_data.cells[L1].pqs[.rho])/dx0
        s1 = (global_data.cells[R0].pqs[.rho] - global_data.cells[L0].pqs[.rho])/dx1
        s2 = (global_data.cells[R1].pqs[.rho] - global_data.cells[R0].pqs[.rho])/dx2
        sL = slope_average(s0, s1)
        sR = slope_average(s1, s2)
        face.left[.rho] = global_data.cells[L0].pqs[.rho] + sL*dxL
        face.right[.rho] = global_data.cells[R0].pqs[.rho] - sR*dxR
        //fmt.printfln("s0= %v s1= %v s2= %v sL= %v sR= %v", s0, s1, s2, sL, sR)
        // p
        s0 = (global_data.cells[L0].pqs[.p] - global_data.cells[L1].pqs[.p])/dx0
        s1 = (global_data.cells[R0].pqs[.p] - global_data.cells[L0].pqs[.p])/dx1
        s2 = (global_data.cells[R1].pqs[.p] - global_data.cells[R0].pqs[.p])/dx2
        sL = slope_average(s0, s1)
        sR = slope_average(s1, s2)
        face.left[.p] = global_data.cells[L0].pqs[.p] + sL*dxL
        face.right[.p] = global_data.cells[R0].pqs[.p] - sR*dxR
        // xvel
        s0 = (global_data.cells[L0].pqs[.xvel] - global_data.cells[L1].pqs[.xvel])/dx0
        s1 = (global_data.cells[R0].pqs[.xvel] - global_data.cells[L0].pqs[.xvel])/dx1
        s2 = (global_data.cells[R1].pqs[.xvel] - global_data.cells[R0].pqs[.xvel])/dx2
        sL = slope_average(s0, s1)
        sR = slope_average(s1, s2)
        face.left[.xvel] = global_data.cells[L0].pqs[.xvel] + sL*dxL
        face.right[.xvel] = global_data.cells[R0].pqs[.xvel] - sR*dxR
        // yvel
        s0 = (global_data.cells[L0].pqs[.yvel] - global_data.cells[L1].pqs[.yvel])/dx0
        s1 = (global_data.cells[R0].pqs[.yvel] - global_data.cells[L0].pqs[.yvel])/dx1
        s2 = (global_data.cells[R1].pqs[.yvel] - global_data.cells[R0].pqs[.yvel])/dx2
        sL = slope_average(s0, s1)
        sR = slope_average(s1, s2)
        face.left[.yvel] = global_data.cells[L0].pqs[.yvel] + sL*dxL
        face.right[.yvel] = global_data.cells[R0].pqs[.yvel] - sR*dxR
        // zvel
        s0 = (global_data.cells[L0].pqs[.zvel] - global_data.cells[L1].pqs[.zvel])/dx0
        s1 = (global_data.cells[R0].pqs[.zvel] - global_data.cells[L0].pqs[.zvel])/dx1
        s2 = (global_data.cells[R1].pqs[.zvel] - global_data.cells[R0].pqs[.zvel])/dx2
        sL = slope_average(s0, s1)
        sR = slope_average(s1, s2)
        face.left[.zvel] = global_data.cells[L0].pqs[.zvel] + sL*dxL
        face.right[.zvel] = global_data.cells[R0].pqs[.zvel] - sR*dxR
        // Fill out other thermo quantities
        face.left[.T] = face.left[.p]/(face.left[.rho]*globals.R_gas)
        face.left[.e] = face.left[.T]*globals.R_gas/(globals.gamma - 1.0)
        face.right[.T] = face.right[.p]/(face.right[.rho]*globals.R_gas)
        face.right[.e] = face.right[.T]*globals.R_gas/(globals.gamma - 1.0)
    }
}


// flux reconstruction via biased averaging procedure
B :: proc (x: complex128) -> complex128 {
    return x/(cmplx.sqrt(1. + x*x) + 1.)
}

B_inv :: proc (x: complex128) -> complex128 {
    return (2.*x)/(1. - x*x)
}

slope_average :: proc (s_l, s_r : complex128) -> complex128 {
    return B_inv((B(s_l) + B(s_r))/2.)
}
