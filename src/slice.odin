package pingvin

import "core:fmt"
import "core:slice"
import "core:os"

Slice :: struct {
    n_cells    :      int,
    first_cell :      Cell_Id,
    last_cell  :      Cell_Id,
    first_x_face:     Interface_Id,
    last_x_face:      Interface_Id,
    interior_faces:   [dynamic]Interface_Id,
    wall_faces:       [dynamic]Interface_Id,
    symm_faces:       [dynamic]Interface_Id,
    up_faces:         #soa[]Interface,
    dn_faces:         #soa[]Interface,
}

Slice_Id :: distinct int

assemble_initial_upstream_interfaces :: proc(up_g: ^Grid_2d) {
    n_faces := len(up_g.quads)
    // upstream faces (i_minus)
    global_data.m_faces = make(#soa[dynamic]Interface, n_faces)
    for i in 0..<n_faces {
        area, normal, t1, t2 := quad_properties(global_data.quads[i])
        global_data.m_faces[i].area = complex(area, 0.0)
        global_data.m_faces[i].normal = normal
        global_data.m_faces[i].t1 = t1
        global_data.m_faces[i].t2 = t2
    }
}

make_face_tag :: proc(q: Quad) -> (tag: string) {
    q := q
    slice.sort(q[:])
    tag = fmt.tprintf("%d-%d-%d-%d", q[0], q[1], q[2], q[3])
    return tag
}

find_face :: proc(needle: Face2, haystack: ^[dynamic]Face2) -> bool {
    for face, i in haystack {
        if faces_are_the_same(needle, face) {
            unordered_remove(haystack, i)
            return true
        }
    }
    return false
}

// It's convenient to assemble cells and slice interfaces at the same
// time because we have that information handy.
assemble_slice_cells_and_interfaces :: proc(slice: ^Slice, up_q, dn_q: []Quad, slice_no: Slice_Id) {
    exists : bool
    jm_id, jp_id, km_id, kp_id : Interface_Id
    n_cells := len(up_q)
    slice.n_cells = n_cells
    slice.first_cell = Cell_Id(int(slice_no)*n_cells)
    slice.last_cell = slice.first_cell + Cell_Id(n_cells) - 1
    slice.first_x_face = Interface_Id(len(global_data.x_faces))
    // We'll build a list of boundary faces and then check off that
    // we've collected them during cell assembly.
    wall_faces : [dynamic]Face2
    append(&wall_faces, ..global_data.wall_boundary.faces[:])
    defer delete(wall_faces)
    symm_faces : [dynamic]Face2
    append(&symm_faces, ..global_data.symm_boundary.faces[:])
    defer delete(symm_faces)
    x_face_tags := make(map[string]Interface_Id)
    defer delete(x_face_tags)
    for i in 0..<n_cells {
        qu := up_q[i]
        qd := dn_q[i]
        hex := hex_from_quads(qu, qd)
        vol := hex_volume(hex)
        ctr := hex_centroid(hex)
        cell_id := Cell_Id(len(global_data.cells))
        cell := Cell{id=cell_id, volume=vol, centroid=ctr}
        append(&global_data.cells, cell)

        // i_minus face (should exist)
        im_id := Interface_Id(int(slice_no)*n_cells + i)
        im_face := global_data.m_faces[im_id]
        global_data.cells[cell_id].faces[.i_minus] = im_id
        global_data.cells[cell_id].outsigns[.i_minus] = complex(1, 0)
        global_data.m_faces[im_id].left_cells = {cell_id, -1}
        if (global_data.m_faces[im_id].right_cells[0] >= 0) {
            // We can set the second left cell of the previous cell.
            prev_cell_id := global_data.m_faces[im_id].right_cells[0]
            prev_face_id := global_data.cells[prev_cell_id].faces[.i_minus]
            global_data.m_faces[prev_face_id].left_cells[1] = cell_id
        }

        // i_plus face (should NOT exist)
        ip_id := Interface_Id((int(slice_no)+1)*n_cells + i)
        area, normal, t1, t2 := quad_properties(qd)
        ip_face := Interface{area=complex(area, 0.0), normal=normal, t1=t1, t2=t2}
        ip_face.right_cells = {cell_id, im_face.right_cells[0]}
        ip_face.left_cells = {-1, -1}
        global_data.cells[cell_id].faces[.i_plus] = ip_id
        global_data.cells[cell_id].outsigns[.i_plus] = complex(-1, 0)
        append(&global_data.m_faces, ip_face)
        
        // j_minus face
        jm_quad := Quad{qu[0], qd[0], qd[1], qu[1]}
        jm_tag := make_face_tag(jm_quad)
        jm_id, exists = x_face_tags[jm_tag]
        if !exists {
            // face does not yet exist
            area, normal, t1, t2 := quad_properties(jm_quad)
            jm_face := Interface{area=complex(area, 0), normal=normal, t1=t1, t2=t2}
            jm_id = Interface_Id(len(global_data.x_faces))
            x_face_tags[jm_tag] = jm_id
            // Set left cell while we know it, and
            // set all other cells to -1 at this point
            // to indicate we need to look later on
            jm_face.left_cells = {cell_id, -1}
            jm_face.right_cells = {-1, -1}
            append(&global_data.x_faces, jm_face)
            global_data.cells[cell_id].outsigns[.j_minus] = complex(1, 0)
            // Since this has just been created, try to figure out
            // what type of interface it is.
            jm_type := InterfaceType.interior
            if slice_no == 0 {
                needle := Face2{qu[0], qu[1]}
                if find_face(needle, &wall_faces) {
                    jm_type = InterfaceType.wall
                }
                else if find_face(needle, &symm_faces) {
                    jm_type = InterfaceType.symm
                }
                #partial switch (jm_type) {
                case .interior:
                    append(&slice.interior_faces, jm_id)
                case .wall:
                    append(&slice.wall_faces, jm_id)
                case .symm:
                    append(&slice.symm_faces, jm_id)
                }
            }
        }
        else {
            // face does exist, so set its right cell
            global_data.x_faces[jm_id].right_cells[0] = cell_id
            global_data.cells[cell_id].outsigns[.j_minus] = complex(-1, 0)
        }
        global_data.cells[cell_id].faces[.j_minus] = jm_id

        // j_plus face
        jp_quad := Quad{qu[2], qd[2], qd[3], qu[3]}
        jp_tag := make_face_tag(jp_quad)
        jp_id, exists = x_face_tags[jp_tag]
        if !exists {
            // face does not yet exist
            area, normal, t1, t2 := quad_properties(jp_quad)
            jp_face := Interface{area=complex(area, 0), normal=normal, t1=t1, t2=t2}
            jp_id = Interface_Id(len(global_data.x_faces))
            x_face_tags[jp_tag] = jp_id
            // Set left cell while we know it, and
            // set all other cells to -1 at this point
            // to indicate we need to look later on
            jp_face.left_cells = {cell_id, -1}
            jp_face.right_cells = {-1, -1}
            append(&global_data.x_faces, jp_face)
            global_data.cells[cell_id].outsigns[.j_plus] = complex(1, 0)
            // Since this has just been created, try to figure out
            // what type of interface it is.
            jp_type := InterfaceType.interior
            if slice_no == 0 {
                needle := Face2{qu[2], qu[3]}
                if find_face(needle, &wall_faces) {
                    jp_type = InterfaceType.wall
                }
                else if find_face(needle, &symm_faces) {
                    jp_type = InterfaceType.symm
                }
                #partial switch (jp_type) {
                case .interior:
                    append(&slice.interior_faces, jp_id)
                case .wall:
                    append(&slice.wall_faces, jp_id)
                case .symm:
                    append(&slice.symm_faces, jp_id)
                }
            }
        }
        else {
            // face does exist, so set its right cell
            global_data.x_faces[jp_id].right_cells[0] = cell_id
            global_data.cells[cell_id].outsigns[.j_plus] = complex(-1, 0)
        }
        global_data.cells[cell_id].faces[.j_plus] = jp_id

        // k_minus face
        km_quad := Quad{qu[3], qd[3], qd[0], qu[0]}
        km_tag := make_face_tag(km_quad)
        km_id, exists = x_face_tags[km_tag]
        if !exists {
            // face does not yet exist
            area, normal, t1, t2 := quad_properties(km_quad)
            km_face := Interface{area=complex(area, 0), normal=normal}
            km_id = Interface_Id(len(global_data.x_faces))
            x_face_tags[km_tag] = km_id
            // Set left cell while we know it, and
            // set all other cells to -1 at this point
            // to indicate we need to look later on
            km_face.left_cells = {cell_id, -1}
            km_face.right_cells = {-1, -1}
            append(&global_data.x_faces, km_face)
            global_data.cells[cell_id].outsigns[.k_minus] = complex(1, 0)
            // Since this has just been created, try to figure out
            // what type of interface it is.
            km_type := InterfaceType.interior
            if slice_no == 0 {
                needle := Face2{qu[3], qu[0]}
                if find_face(needle, &wall_faces) {
                    km_type = InterfaceType.wall
                }
                else if find_face(needle, &symm_faces) {
                    km_type = InterfaceType.symm
                }
                #partial switch (km_type) {
                case .interior:
                    append(&slice.interior_faces, km_id)
                case .wall:
                    append(&slice.wall_faces, km_id)
                case .symm:
                    append(&slice.symm_faces, km_id)
                }
            }
        }
        else {
            // face does exist, so set its right cell
            global_data.x_faces[km_id].right_cells[0] = cell_id
            global_data.cells[cell_id].outsigns[.k_minus] = complex(-1, 0)
        }
        global_data.cells[cell_id].faces[.k_minus] = km_id

        // k_plus face
        kp_quad := Quad{qu[1], qd[1], qd[2], qu[2]}
        kp_tag := make_face_tag(kp_quad)
        kp_id, exists = x_face_tags[kp_tag]
        if !exists {
            // face does not yet exist
            area, normal, t1, t2 := quad_properties(kp_quad)
            kp_face := Interface{area=complex(area, 0), normal=normal, t1=t1, t2=t2}
            kp_id = Interface_Id(len(global_data.x_faces))
            x_face_tags[kp_tag] = kp_id
            // Set left cell while we know it, and
            // set all other cells to -1 at this point
            // to indicate we need to look later on
            kp_face.left_cells = {cell_id, -1}
            kp_face.right_cells = {-1, -1}
            append(&global_data.x_faces, kp_face)
            global_data.cells[cell_id].outsigns[.k_plus] = complex(1, 0)
            // Since this has just been created, try to figure out
            // what type of interface it is.
            kp_type := InterfaceType.interior
            if slice_no == 0 {
                needle := Face2{qu[1], qu[2]}
                if find_face(needle, &wall_faces) {
                    kp_type = InterfaceType.wall
                }
                else if find_face(needle, &symm_faces) {
                    kp_type = InterfaceType.symm
                }
                #partial switch (kp_type) {
                case .interior:
                    append(&slice.interior_faces, kp_id)
                case .wall:
                    append(&slice.wall_faces, kp_id)
                case .symm:
                    append(&slice.symm_faces, kp_id)
                }
            }
        }
        else {
            // face does exist, so set its right cell
            global_data.x_faces[kp_id].right_cells[0] = cell_id
            global_data.cells[cell_id].outsigns[.k_plus] = complex(-1, 0)
        }
        global_data.cells[cell_id].faces[.k_plus] = kp_id
    }
    slice.last_x_face = Interface_Id(len(global_data.x_faces)-1)
    // With all cells and interfaces assembled, we can build out
    // the stencils for the interfaces in the cross-stream direction
    for i in slice.first_x_face..=slice.last_x_face {
        // Fill stencil to left
        lc_id := global_data.x_faces[i].left_cells[0]
        if lc_id >= 0 {
            // find this interface in that cell's list...
            switch (i) {
            case global_data.cells[lc_id].faces[.j_minus]:
                id := global_data.cells[lc_id].faces[.j_plus]
                global_data.x_faces[i].left_cells[1] = global_data.x_faces[id].left_cells[0]
            case global_data.cells[lc_id].faces[.j_plus]:
                id := global_data.cells[lc_id].faces[.j_minus]
                global_data.x_faces[i].left_cells[1] = global_data.x_faces[id].left_cells[0]
            case global_data.cells[lc_id].faces[.k_minus]:
                id := global_data.cells[lc_id].faces[.k_plus]
                global_data.x_faces[i].left_cells[1] = global_data.x_faces[id].left_cells[0]
            case global_data.cells[lc_id].faces[.k_plus]:
                id := global_data.cells[lc_id].faces[.k_minus]
                global_data.x_faces[i].left_cells[1] = global_data.x_faces[id].left_cells[0]
            }
        }
        // Fill stencil to right
        rc_id := global_data.x_faces[i].right_cells[0]
        if rc_id >= 0 {
            // find this interface in that cell's list...
            switch (i) {
            case global_data.cells[rc_id].faces[.j_minus]:
                id := global_data.cells[rc_id].faces[.j_plus]
                global_data.x_faces[i].right_cells[1] = global_data.x_faces[id].right_cells[0]
            case global_data.cells[rc_id].faces[.j_plus]:
                id := global_data.cells[rc_id].faces[.j_minus]
                global_data.x_faces[i].right_cells[1] = global_data.x_faces[id].right_cells[0]
            case global_data.cells[rc_id].faces[.k_minus]:
                id := global_data.cells[rc_id].faces[.k_plus]
                global_data.x_faces[i].right_cells[1] = global_data.x_faces[id].right_cells[0]
            case global_data.cells[rc_id].faces[.k_plus]:
                id := global_data.cells[rc_id].faces[.k_minus]
                global_data.x_faces[i].right_cells[1] = global_data.x_faces[id].right_cells[0]
            }
        }
    }
    // When building initial slice, we should now check that all boundary faces were assigned.
    if slice_no == 0 {
        if len(wall_faces) != 0 {
            fmt.printfln("ERROR: There are unassigned wall faces. Number of unassigned faces are: %v", len(wall_faces))
            fmt.println("Those faces are:")
            for f in wall_faces {
                fmt.println(f)
            }
            fmt.println("Exiting!")
            os.exit(1)
        }
        if len(symm_faces) != 0 {
            fmt.printfln("ERROR: There are unassigned symm faces. Number of unassigned faces are: %v", len(symm_faces))
            fmt.println("Those faces are:")
            for f in symm_faces {
                fmt.println(f)
            }
            fmt.println("Exiting!")
            os.exit(1)
        }
    }
    else {
        // For all subsequent slices, we can easily construct the arrays of handles by incrementing indices
        slice0 := global_data.slices[0]
        n_faces := len(slice0.interior_faces) + len(slice0.wall_faces) + len(slice0.symm_faces)
        append(&slice.interior_faces, ..slice0.interior_faces[:])
        for &face in slice.interior_faces {
            face += Interface_Id(n_faces*int(slice_no))
        }
        append(&slice.wall_faces, ..slice0.wall_faces[:])
        for &face in slice.wall_faces {
            face += Interface_Id(n_faces*int(slice_no))
        }
        append(&slice.symm_faces, ..slice0.symm_faces[:])
        for &face in slice.symm_faces {
            face += Interface_Id(n_faces*int(slice_no))
        }
    }
    // Set slices to upstream and downstream interfaces
    slice.up_faces = global_data.m_faces[slice.first_cell:slice.last_cell+1]
    slice.dn_faces = global_data.m_faces[slice.last_cell:slice.last_cell+Cell_Id(n_cells)+1]
}

create_slice :: proc (x: f64, xsect: ^Cross_Section, slice: Slice_Id) {
    up_g := &global_data.up_grid
    dn_g := &global_data.dn_grid
    compute_grid_2d(dn_g, up_g, &global_data.rtheta_grid, xsect)
    add_3d_slice_of_hexes(up_g, dn_g)
    if (slice == 0) {
        assemble_initial_upstream_interfaces(up_g)
    }
    append(&global_data.slices, Slice{})
    assemble_slice_cells_and_interfaces(&global_data.slices[slice], up_g.quads[:], dn_g.quads[:], slice)
    if (slice == 0) {
        apply_inflow(&global_data.slices[0])
    }
    else {
        prep_slice(&global_data.slices[slice], &global_data.slices[slice-1])
    }
}

update_primitives :: proc (slice: ^Slice) {
    for i in slice.first_cell..=slice.last_cell {
        global_data.cells[i].pqs = prim_from_cq(global_data.cells[i].cqs)
    }
}

compute_interior_fluxes :: proc (slice: ^Slice) {
    for f_id in slice.interior_faces {
        f := &global_data.x_faces[f_id]
        f.flux = flux_calc(f.left, f.right)
    }
}

compute_residuals :: proc (slice: ^Slice) {
    for &cell in global_data.cells[slice.first_cell:slice.last_cell+1] {
        compute_residual(&cell)
    }
}

eval_slice_residual :: proc (slice: ^Slice) {
    // 1. Update primitive values
    update_primitives(slice)

    // 2. Apply reconstruction (at interior and boundaries)
    low_order_reconstruction(slice.interior_faces[:])
    low_order_recon_boundary(slice.wall_faces[:])
    low_order_recon_boundary(slice.symm_faces[:])
    low_order_recon_downstream(slice.dn_faces)

    // 3. Compute fluxes
    compute_interior_fluxes(slice)
    apply_slip_wall_flux(slice.wall_faces[:])
    apply_symm_flux(slice.symm_faces[:])
    apply_downstream_flux(slice.dn_faces)

    // 4. Compute residual per cell
    compute_residuals(slice)
}

prep_slice :: proc (slice, prev_slice: ^Slice) {
    curr_offset := slice.first_cell
    prev_offset := prev_slice.first_cell
    for i in 0..<slice.n_cells {
        i := Cell_Id(i)
        global_data.cells[i+curr_offset].cqs = global_data.cells[i+prev_offset].cqs
        global_data.cells[i+curr_offset].pqs = global_data.cells[i+prev_offset].pqs
    }
}


