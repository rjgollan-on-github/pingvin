package pingvin

import "core:math"
import "core:fmt"
import "core:os"

ZERO_TOL :: 1.0e-12

GMRES_Workspace :: struct {
    nvars:             int,
    // Vectors with length of solution vector
    r0, rhs, v, w, x:  []f64,
    // VT is Krylov subspace stacked into flat array
    VT:                []f64,
    // Storage for GMRES support
    h, hR:             []f64,
    g0, g1:            []f64,
    H0, H1:            [][]f64,
    Q0, Q1:            [][]f64,
    Omega:             [][]f64,
}

gws := GMRES_Workspace{}

allocate_GMRES_Workspace :: proc (nvars, max_iterations : int) {
    gws.nvars = nvars
    gws.r0 = make([]f64, nvars)
    gws.rhs = make([]f64, nvars)
    gws.v = make([]f64, nvars)
    gws.w = make([]f64, nvars)
    gws.x = make([]f64, nvars)
    gws.VT = make([]f64, nvars*(max_iterations+1))
    gws.h = make([]f64, max_iterations+1)
    gws.hR = make([]f64, max_iterations+1)
    gws.g0 = make([]f64, max_iterations+1)
    gws.g1 = make([]f64, max_iterations+1)
    gws.H0 = make([][]f64, max_iterations+1)
    for i in 0..<max_iterations+1 do gws.H0[i] = make([]f64, max_iterations)
    gws.H1 = make([][]f64, max_iterations+1)
    for i in 0..<max_iterations+1 do gws.H1[i] = make([]f64, max_iterations)
    gws.Q0 = make([][]f64, max_iterations+1)
    for i in 0..<max_iterations+1 do gws.Q0[i] = make([]f64, max_iterations+1)
    gws.Q1 = make([][]f64, max_iterations+1)
    for i in 0..<max_iterations+1 do gws.Q1[i] = make([]f64, max_iterations+1)
    gws.Omega = make([][]f64, max_iterations+1)
    for i in 0..<max_iterations+1 do gws.Omega[i] = make([]f64, max_iterations+1)
    return
}

delete_GMRES_Workspace :: proc() {
    delete(gws.r0)
    delete(gws.rhs)
    delete(gws.v)
    delete(gws.w)
    delete(gws.x)
    delete(gws.VT)
    delete(gws.h)
    delete(gws.hR)
    delete(gws.g0)
    delete(gws.g1)
    for i in 0..<len(gws.H0) do delete(gws.H0[i])
    delete(gws.H0)
    for i in 0..<len(gws.H1) do delete(gws.H1[i])
    delete(gws.H1)
    for i in 0..<len(gws.Q0) do delete(gws.Q0[i])
    delete(gws.Q0)
    for i in 0..<len(gws.Q1) do delete(gws.Q1[i])
    delete(gws.Q1)
    for i in 0..<len(gws.Omega) do delete(gws.Omega[i])
    delete(gws.Omega)
    return
}

slice_residual_norm :: proc(slice: ^Slice) -> f64 {
    l2norm := 0.0
    for i in slice.first_cell..=slice.last_cell {
        cell := global_data.cells[i]
        for cq in cell.R {
            l2norm += real(cq)*real(cq)
        }
    }
    l2norm = math.sqrt(l2norm)
    return l2norm
}

slice_dU_norm :: proc(slice: ^Slice) -> f64 {
    l2norm := 0.0
    for i in slice.first_cell..=slice.last_cell {
        cell := global_data.cells[i]
        for cq in cell.dU {
            l2norm += cq*cq
        }
    }
    l2norm = math.sqrt(l2norm)
    return l2norm
}

update_slice :: proc(slice: ^Slice) {
    for i in slice.first_cell..=slice.last_cell {
        cell := &global_data.cells[i]
 //       fmt.println("Before update")
 //       fmt.println(cell.cqs0)
        cell.cqs0 += cell.dU
//        fmt.println("After update")
//        fmt.println(cell.cqs0)
        for cq in Conserved_Quantities {
            cell.cqs[cq] = complex(cell.cqs0[cq], 0.0)
        }
    }
}

eval_jacobian_vector_product :: proc (slice: ^Slice, v, Jv: []f64) {
    sigma := globals.cfg.perturbation_size
    n_cons := len(Conserved_Quantities)
    // Place perturbed quantity in cqs
    for &cell, i in global_data.cells[slice.first_cell:slice.last_cell+1] {
        for cq, icq in Conserved_Quantities {
            cell.cqs[cq] = complex(cell.cqs0[cq], 0) + complex(0, sigma*v[i*n_cons + icq])
        }
    }
    // Evaluate residual
    eval_slice_residual(slice)
    // Compute Jv using Frechet derivative (in complex space)
    for &cell, i in global_data.cells[slice.first_cell:slice.last_cell+1] {
        for cq, icq in Conserved_Quantities {
            Jv[i*n_cons + icq] = imag(cell.R[cq])/sigma
        }
    }
    return
}

solve_slice :: proc (slice_no: Slice_Id) -> (is_converged : bool) {
    slice := &global_data.slices[slice_no]
    cfg := globals.cfg
    print_every_n_slice := cfg.print_every_n_slice
    is_converged = false
    // 0. Determine residual norm at start
    eval_slice_residual(slice)
    R0_norm := slice_residual_norm(slice)
    copy_slice_to_base_state(slice)
    copy_slice_residual_to_gws(slice)

    if R0_norm < cfg.slice_absolute_residual {
        is_converged = true
        return is_converged
    }

    // 1. Iterations: continue until a desired convergence, or max Newton steps
    max_steps := cfg.max_newton_steps
    for step in 0..<max_steps {
        // Solve
        is_gmres_converged := gmres_solve(slice)
        // Update
        update_slice(slice)
        eval_slice_residual(slice)
        R_norm := slice_residual_norm(slice)
        copy_slice_to_base_state(slice)
        copy_slice_residual_to_gws(slice)
        // Check on convergence
        rel_residual := R_norm/R0_norm
        dU_norm := slice_dU_norm(slice)
        if (int(slice_no) % print_every_n_slice) == 0 {
            fmt.printfln("     [slice-%03d]: step= %d, rel. residual= %.6e  ||dU||= %.6e", slice_no, step, rel_residual, dU_norm)
        }
        if rel_residual < cfg.slice_relative_residual {
            is_converged = true
            if (int(slice_no) % print_every_n_slice) == 0 {
                fmt.println("       !!! Slice converged: relative residual target achieved. !!!")
            }
            return is_converged
        }
        if dU_norm < cfg.slice_change_in_update {
            is_converged = true
            if (int(slice_no) % print_every_n_slice) == 0 {
                fmt.println("       !!! Slice converged: target for change over step (||dU||) achieved. !!!")
            }
            return is_converged
        }
    }
    return is_converged
}

copy_slice_to_base_state :: proc (slice: ^Slice) {
    for &cell, i in global_data.cells[slice.first_cell:slice.last_cell+1] {
        copy_residual_to_base_state(&cell);
        copy_cqs_to_base_state(&cell);
    }
}

copy_slice_residual_to_gws :: proc (slice: ^Slice) {
    n_cons := len(Conserved_Quantities)
    for cell, i in global_data.cells[slice.first_cell:slice.last_cell+1] {
        for cq, icq in Conserved_Quantities {
            gws.rhs[i*n_cons + icq] = -cell.Ru[cq]
        }
    }
    return
}

compute_r0 :: proc () {
    // Presently, we assume x0[] = 0.0, so r0 = rhs
    copy(gws.r0, gws.rhs)
}

copy_matrix :: proc (tgt, src: [][]f64) {
    for i in 0..<len(tgt) {
        for j in 0..<len(tgt[i]) do tgt[i][j] = src[i][j]
    }
}

copy_array :: proc (tgt, src: []f64) {
    for i in 0..<len(tgt) do tgt[i] = src[i]
    return
}

l2norm :: proc (vec: []f64) -> f64 {
    sum : f64 = 0.0
    for v in vec {
        sum += v*v
    }
    return math.sqrt(sum)
}

dot_array_f64 :: proc(a, b: []f64) -> f64 {
    sum : f64  = 0.0
    for i in 0..<len(a) {
        sum += a[i]*b[i]
    }
    return sum
}

dot_matrix_vector :: proc(A: [][]f64, a_row, a_col: int, b, c: []f64) {
    for i in 0..<len(c) do c[i] = 0.0
    for i in 0..=a_row {
        for j in 0..=a_col {
            c[i] += A[i][j] * b[j]
        }
    }
    return
}

dot_matrix_upper :: proc(A: [][]f64, a_row, a_col: int, B: [][]f64, b_col: int, C: [][]f64) {
    for r in 0..=a_row {
        for c in 0..=b_col {
            C[r][c] = 0.0
            for i in 0..=a_col {
                C[r][c] += A[r][i] * B[i][c]
            }
        }
    }
    return
}

// Collect the dot operations in a group
dot :: proc {
    dot_vector3,
    dot_array_f64,
    dot_matrix_vector,
    dot_matrix_upper,
}

upper_solve :: proc (A: [][]f64, n: int, b : []f64) {
    b[n-1] /= A[n-1][n-1]
    for i := n-2; i >= 0; i -= 1 {
        sum := b[i]
        for j in i+1..<n do sum -= A[i][j] * b[j]
        b[i] = sum/A[i][i]
    }
    return
}

transpose_and_dot :: proc (A: []f64, stride, a_row, a_col : int, b, c: []f64) {
    for i in 0..<len(c) do c[i] = 0.0
    for i in 0..=a_row {
        for j in 0..=a_col {
            c[j] += A[i*stride+j] * b[i]
        }
    }
}


eye :: proc (A: [][]f64) {
    for i in 0..<len(A) {
        for j in 0..<len(A[i]) {
            A[i][j] = 0.0
        }
        A[i][i] = 1.0
    }
    return
}

zeros_matrix :: proc (A: [][]f64) {
    for i in 0..<len(A) {
        for j in 0..<len(A[i]) {
            A[i][j] = 0.0
        }
    }
    return
}

zeros_vector :: proc (a: []f64) {
    for i in 0..<len(a) do a[i] = 0.0
}

zeros :: proc {
    zeros_matrix,
    zeros_vector,
}

gmres_solve :: proc (slice: ^Slice) -> (is_converged : bool) {
    n_cons := len(Conserved_Quantities)
    m, ok := perform_gmres_iterations(slice)
    is_converged = ok

    upper_solve(gws.H1, m, gws.g1)

    transpose_and_dot(gws.VT, gws.nvars, m-1, gws.nvars-1, gws.g1, gws.x)
    for &cell, i in global_data.cells[slice.first_cell:slice.last_cell+1] {
        for cq, icq in Conserved_Quantities {
            cell.dU[cq] = gws.x[i*n_cons + icq]
        }
    }
    return is_converged
}

multiply :: proc (A: [][]f64, x, b : []f64) {
    for i in 0..<len(A) {
        b[i] = 0.0
        for j in 0..<len(A[i]) {
            b[i] += A[i][j]*x[j]
        }
    }
    return
}

perform_gmres_iterations :: proc (slice: ^Slice) -> (n_iterations: int, is_converged : bool) {
    cfg := globals.cfg
    is_converged = false
    nvars := gws.nvars
    m := cfg.max_gmres_iterations
    // 0. Set support matrices and vectors for GMRES
    zeros(gws.H0)
    zeros(gws.H1)
    eye(gws.Omega)
    zeros(gws.g0)
    zeros(gws.g1)
    // 1. Compute r0 = b - Ax0; beta := ||r0||2, and v1 := r0/beta
    compute_r0();
    beta0 := l2norm(gws.r0) // store for diagnostics
    beta := beta0
    gws.g0[0] = beta0
    for i in 0..<nvars {
        gws.v[i] = gws.r0[i]/beta
        gws.VT[i] = gws.v[i]
    }
    target_residual := cfg.gmres_relative_residual * beta

    //fmt.printfln("gmres: beta= %.16e, r0[0]= %.16e", beta, gws.r0[0])
    //fmt.printfln("gmres: target_residual= %.16e", target_residual)

    // 2. begin loop
    for j in 0..<cfg.max_gmres_iterations {
//    for j in 0..<1 {
        // 3. Compute w := Jv
        eval_jacobian_vector_product(slice, gws.v, gws.w)
        //multiply(A, gws.v, gws.w)

        // 4. Begin orthonormalisation
        for i in 0..=j {
            // 5. h_ij := (wj,vi)
            // Extract v from VT
            offset := i*nvars
            for ii in 0..<nvars do gws.v[ii] = gws.VT[ii + offset]
            // Now ready for dot product
            h_ij := dot(gws.w, gws.v)
            gws.H0[i][j] = h_ij
            // 6. w_j := w_j - h_ij*v_i
            for ii in 0..<nvars do gws.w[ii] -= h_ij*gws.v[ii]

        } // 7. end 4.
        // 8. 
        h_jp1j := l2norm(gws.w)
        gws.H0[j+1][j] = h_jp1j
//        fmt.printfln("gmres: h0[%d,%d]= %v", j+1, j, h_jp1j)
        // 9.
        offset := (j+1)*nvars
        for ii in 0..<nvars {
            gws.v[ii] = gws.w[ii]/h_jp1j
            gws.VT[ii+offset] = gws.v[ii]
        }

        // Build rotated Hessenberg matrix progressively
        if j != 0 {
            // Extract final column in H
            for i in 0..=j do gws.h[i] = gws.H0[i][j]
            // Rotate column by previous rotations (stored in Q)
            dot(gws.Q0, j, j, gws.h, gws.hR)
            // Place column back in H
            for i in 0..=j do gws.H0[i][j] = gws.hR[i]
        }

        // Form new Omega
        eye(gws.Omega)
        c_j, s_j, denom : f64
        denom = math.sqrt(gws.H0[j][j]*gws.H0[j][j] + gws.H0[j+1][j]*gws.H0[j+1][j])
        s_j = gws.H0[j+1][j]/denom
        c_j = gws.H0[j][j]/denom
        gws.Omega[j][j] = c_j
        gws.Omega[j][j+1] = s_j
        gws.Omega[j+1][j] = -s_j
        gws.Omega[j+1][j+1] = c_j
        // Apply rotations
        dot_matrix_upper(gws.Omega, j+1, j+1, gws.H0, j, gws.H1)
        dot_matrix_vector(gws.Omega, j+1, j+1, gws.g0, gws.g1)

        // Accumulate rotations in Q
        if j == 0 {
            copy_matrix(gws.Q1, gws.Omega)
        }
        else {
            dot_matrix_upper(gws.Omega, j+1, j+1, gws.Q0, j+1, gws.Q1)
        }
        // Prepare for next step
        copy_matrix(gws.H0, gws.H1)
        copy_matrix(gws.Q0, gws.Q1)
        copy_array(gws.g0, gws.g1)
        
        // Get residual
        residual := abs(gws.g1[j+1])
        if residual <= target_residual {
            is_converged = true
            return j+1, is_converged
        }
    }
    return m, is_converged
}
