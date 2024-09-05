package pingvin

import "core:fmt"
// For interior faces
low_order_reconstruction :: proc (faces: []Interface_Id) {
    for f_id in faces {
        face := &global_data.x_faces[f_id]
        face.left = global_data.cells[face.left_cells[0]].pqs
        face.right = global_data.cells[face.right_cells[0]].pqs
        if f_id == 0 || f_id == 1 || f_id == 2 || f_id == 3 {
            pos := Vector3{}
            for i in face.quad {
                pos += global_data.vertices[i]
            }
            pos *= 0.25
            sigma := 1.0e-250
            fmt.printfln("--- f_id= %d, pos= %v ---", f_id, pos)
            fmt.printfln("L: u= %.8e u-d= %.8e", real(face.left[.xvel]), imag(face.left[.xvel])/sigma)
            fmt.printfln("L: v= %.8e v-d= %.8e", real(face.left[.yvel]), imag(face.left[.yvel])/sigma)
            fmt.printfln("L: w= %.8e w-d= %.8e", real(face.left[.zvel]), imag(face.left[.zvel])/sigma)
            fmt.printfln("L: p= %.8e p-d= %.8e", real(face.left[.p]), imag(face.left[.p])/sigma)
            fmt.printfln("L: rho= %.8e rho-d= %.8e", real(face.left[.rho]), imag(face.left[.rho])/sigma)
            fmt.printfln("R: u= %.8e u-d= %.8e", real(face.right[.xvel]), imag(face.right[.xvel])/sigma)
            fmt.printfln("R: v= %.8e v-d= %.8e", real(face.right[.yvel]), imag(face.right[.yvel])/sigma)
            fmt.printfln("R: w= %.8e w-d= %.8e", real(face.right[.zvel]), imag(face.right[.zvel])/sigma)
            fmt.printfln("R: p= %.8e p-d= %.8e", real(face.right[.p]), imag(face.right[.p])/sigma)
            fmt.printfln("R: rho= %.8e rho-d= %.8e", real(face.right[.rho]), imag(face.right[.rho])/sigma)
        }
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


