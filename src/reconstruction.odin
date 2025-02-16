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
        face.left = global_data.cells[face.left_cells[0]].pqs
    }
}

high_order_recon_downstream :: proc (faces: #soa[]Interface) {
    pqs_for_reconstruction := [?]Primitive_Quantities{ .rho, .p, .xvel, .yvel, .zvel }
    sR := complex(0.0, 0.0)
    for &face in faces {
        L1 := face.left_cells[1]
        L0 := face.left_cells[0]
        dx := magnitude(global_data.cells[L0].centroid - global_data.cells[L1].centroid)
        dxL := magnitude(face.centroid - global_data.cells[L0].centroid)
        for pq in pqs_for_reconstruction {
            sL := (global_data.cells[L0].pqs[pq] - global_data.cells[L1].pqs[pq])/dx
            //s := slope_average(sL, sR)
            face.left[pq] = global_data.cells[L0].pqs[pq] + sL*dxL
        }
        // Fill out other thermo quantities
        face.left[.T] = face.left[.p]/(face.left[.rho]*globals.R_gas)
        face.left[.e] = face.left[.T]*globals.R_gas/(globals.gamma - 1.0)
    }
}

// At slice boundaries
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

high_order_recon_boundary :: proc (faces: []Interface_Id) {
    for f_id in faces {
        face := &global_data.x_faces[f_id]
        if face.left_cells[0] >= 0 {
            L1 := face.left_cells[1]
            L0 := face.left_cells[0]
            dx := magnitude(global_data.cells[L0].centroid - global_data.cells[L1].centroid)
            // Only need to handle pressure
            sL := (global_data.cells[L0].pqs[.p] - global_data.cells[L1].pqs[.p])/dx
            sR := complex(0.0, 0.0)
            s := slope_average(sL, sR)
            dxL := magnitude(face.centroid - global_data.cells[L0].centroid)
            face.left[.p] = global_data.cells[L0].pqs[.p] + s*dxL
        }
        else { // must be domain cell on right
            R1 := face.right_cells[1]
            R0 := face.right_cells[0]
            dx := magnitude(global_data.cells[R1].centroid - global_data.cells[R0].centroid)
            // Only need to handle pressure
            sL := complex(0.0, 0.0)
            sR := (global_data.cells[R1].pqs[.p] - global_data.cells[R0].pqs[.p])/dx
            s := slope_average(sL, sR)
            dxR := magnitude(face.centroid - global_data.cells[R0].centroid)
            face.right[.p] = global_data.cells[R0].pqs[.p] - s*dxR
        }
    }
}

high_order_recon_interior :: proc (faces: []Interface_Id) {
    s0, s1, s2 : complex128
    sL, sR : complex128
    pqs_for_reconstruction := [?]Primitive_Quantities{ .rho, .p, .xvel, .yvel, .zvel }

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

        for pq in pqs_for_reconstruction {
            s0 = (global_data.cells[L0].pqs[pq] - global_data.cells[L1].pqs[pq])/dx0
            s1 = (global_data.cells[R0].pqs[pq] - global_data.cells[L0].pqs[pq])/dx1
            s2 = (global_data.cells[R1].pqs[pq] - global_data.cells[R0].pqs[pq])/dx2
            sL = slope_average(s0, s1)
            sR = slope_average(s1, s2)
            face.left[pq] = global_data.cells[L0].pqs[pq] + sL*dxL
            face.right[pq] = global_data.cells[R0].pqs[pq] - sR*dxR
        }
        // Fill out other thermo quantities
        face.left[.T] = face.left[.p]/(face.left[.rho]*globals.R_gas)
        face.left[.e] = face.left[.T]*globals.R_gas/(globals.gamma - 1.0)
        face.right[.T] = face.right[.p]/(face.right[.rho]*globals.R_gas)
        face.right[.e] = face.right[.T]*globals.R_gas/(globals.gamma - 1.0)
    }
}

right_near_wall_reconstruction :: proc (faces: []Interface_Id) {
    s0, s1, sL : complex128
    pqs_for_reconstruction := [?]Primitive_Quantities{ .rho, .p, .xvel, .yvel, .zvel }

    for f_id in faces {
        face := &global_data.x_faces[f_id]
        L1 := face.left_cells[1]
        L0 := face.left_cells[0]
        R0 := face.right_cells[0]
        dx0 := magnitude(global_data.cells[L0].centroid - global_data.cells[L1].centroid)
        dx1 := magnitude(global_data.cells[R0].centroid - global_data.cells[L0].centroid)
        dxL := magnitude(face.centroid - global_data.cells[L0].centroid)

        for pq in pqs_for_reconstruction {
            s0 = (global_data.cells[L0].pqs[pq] - global_data.cells[L1].pqs[pq])/dx0
            s1 = (global_data.cells[R0].pqs[pq] - global_data.cells[L0].pqs[pq])/dx1
            sL = slope_average(s0, s1)
            face.left[pq] = global_data.cells[L0].pqs[pq] + sL*dxL
            face.right[pq] = global_data.cells[R0].pqs[pq]
        }
        face.left[.T] = face.left[.p]/(face.left[.rho]*globals.R_gas)
        face.left[.e] = face.left[.T]*globals.R_gas/(globals.gamma - 1.0)
        face.right[.T] = face.right[.p]/(face.right[.rho]*globals.R_gas)
        face.right[.e] = face.right[.T]*globals.R_gas/(globals.gamma - 1.0)
    }
}

left_near_wall_reconstruction :: proc (faces: []Interface_Id) {
    s1, s2, sR : complex128
    pqs_for_reconstruction := [?]Primitive_Quantities{ .rho, .p, .xvel, .yvel, .zvel }

    for f_id in faces {
        face := &global_data.x_faces[f_id]
        L0 := face.left_cells[0]
        R0 := face.right_cells[0]
        R1 := face.right_cells[1]
        dx1 := magnitude(global_data.cells[R0].centroid - global_data.cells[L0].centroid)
        dx2 := magnitude(global_data.cells[R1].centroid - global_data.cells[R0].centroid)
        dxR := magnitude(face.centroid - global_data.cells[R0].centroid)

        for pq in pqs_for_reconstruction {
            s1 = (global_data.cells[R0].pqs[pq] - global_data.cells[L0].pqs[pq])/dx1
            s2 = (global_data.cells[R1].pqs[pq] - global_data.cells[R0].pqs[pq])/dx2
            sR = slope_average(s1, s2)
            face.left[pq] = global_data.cells[L0].pqs[pq]
            face.right[pq] = global_data.cells[R0].pqs[pq] - sR*dxR
        }
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
