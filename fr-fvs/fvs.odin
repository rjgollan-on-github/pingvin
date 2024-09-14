package fvs

import "core:math"
import "core:os"
import "core:fmt"

GAMMA :: 1.4
R :: 287.0
FILENAME :: "sod-second-order.dat"

L := 1.0
dx := 0.01
dt := 1.0e-6
nsteps := 600

Conserved_Quantities :: enum {mass, mom, energy}
Primitive_Quantities :: enum {rho, p, u, e}

vl_flux :: proc (P: [Primitive_Quantities]f64) -> (f_plus, f_minus: [Conserved_Quantities]f64) {
	c := math.sqrt(GAMMA*P[.p]/P[.rho])
	M := P[.u]/c
	E := P[.rho]*(P[.e] + 0.5*P[.u]*P[.u])

	if M >= 1.0 {
		f_plus[.mass] = P[.rho]*P[.u]
		f_plus[.mom] = P[.rho]*P[.u]*P[.u] + P[.p]
		f_plus[.energy] = P[.u]*(E + P[.p])
		f_minus[.mass] = 0.0
		f_minus[.mom] = 0.0
		f_minus[.energy] = 0.0
	}
	else if M <= -1.0 {
		f_plus[.mass] = 0.0
		f_plus[.mom] = 0.0
		f_plus[.energy] = 0.0
		f_minus[.mass] = P[.rho]*P[.u]
		f_minus[.mom] = P[.rho]*P[.u]*P[.u] + P[.p]
		f_minus[.energy] = P[.u]*(E + P[.p])
	}
	else {
		f_plus[.mass] = P[.rho]*c*(M + 1.0)*(M + 1.0)/4.0
		f_plus[.mom] = f_plus[.mass]*((GAMMA - 1.0)*P[.u] + 2.0*c)/GAMMA
		f_plus[.energy] = f_plus[.mass]*(math.pow((GAMMA - 1.0)*P[.u] + 2.0*c, 2))/(2.0*(GAMMA*GAMMA - 1.0))
		f_minus[.mass] = -P[.rho]*c*(M - 1.0)*(M - 1.0)/4.0
		f_minus[.mom] = f_minus[.mass]*((GAMMA - 1.0)*P[.u] - 2.0*c)/GAMMA
		f_minus[.energy] = f_minus[.mass]*(math.pow((GAMMA - 1.0)*P[.u] - 2.0*c, 2))/(2.0*(GAMMA*GAMMA - 1.0))
	}
	return f_plus, f_minus
}

/*
 Alternative for bias-averaging function
B :: proc (x: f64) -> f64 {
	return math.atan(x)
}

B_inv :: proc (x: f64) -> f64 {
	return math.tan(x)
}
*/

B :: proc (x: f64) -> f64 {
	return x/(math.sqrt(1.0 + x*x) + 1.0)
}

B_inv :: proc (x: f64) -> f64 {
	return (2.0*x)/(1.0 - x*x)
}


flux_slope :: proc (s_l, s_r : f64) -> f64 {
	return B_inv((B(s_l) + B(s_r))/2.0)
}

main :: proc () {
	ncells := int(L/dx)
	// Set up x locations
	xs : []f64
	xs = make([]f64, ncells)
	defer delete(xs)
	for i in 0..<ncells {
		xs[i] = (f64(i)+0.5)*dx
	}
	// Set up arrays for conserved quantities and primitive
	U : [][Conserved_Quantities]f64
	U = make([][Conserved_Quantities]f64, ncells)
	defer delete(U)
	
	P : [][Primitive_Quantities]f64
	P = make([][Primitive_Quantities]f64, ncells)
	defer delete(P)
	
	f_p : [][Conserved_Quantities]f64
	f_p = make([][Conserved_Quantities]f64, ncells)
	defer delete(f_p)
	
	f_m : [][Conserved_Quantities]f64
	f_m = make([][Conserved_Quantities]f64, ncells)
	defer delete(f_m)

	s_p : [][Conserved_Quantities]f64
	s_p = make([][Conserved_Quantities]f64, ncells)
	defer delete(s_p)

	s_m : [][Conserved_Quantities]f64
	s_m = make([][Conserved_Quantities]f64, ncells)
	defer delete(s_m)
		
	for i in 0..<ncells {
		x := xs[i]
		if x < 0.5 {
	      // high pressure condition
    	  P[i][.p] = 1.0e5
	      P[i][.rho] = 1.0
    	}
		else {
	      P[i][.p] = 1.0e4
    	  P[i][.rho] = 0.125
		}
	    P[i][.u] = 0.0
	    P[i][.e] = P[i][.p]/((GAMMA - 1.0)*P[i][.rho])
	    // conserved quantities
	    U[i][.mass] = P[i][.rho]
	    U[i][.mom] = P[i][.rho]*P[i][.u]
	    U[i][.energy] = P[i][.rho]*(P[i][.e] + 0.5*P[i][.rho]*P[i][.u]*P[i][.u])
    }

	F : [][Conserved_Quantities]f64
	F = make([][Conserved_Quantities]f64, ncells+1)
	// Set end fluxes
	F[0][.mass] = P[0][.rho]*P[0][.u]
	F[0][.mom] = P[0][.rho]*P[0][.u]*P[0][.u] + P[0][.p]
	F[0][.energy] = P[0][.u]*(P[0][.rho]*(P[0][.e] + 0.5*P[0][.u]*P[0][.u]) + P[0][.p])
	F[ncells][.mass] = P[ncells-1][.rho]*P[ncells-1][.u]
	F[ncells][.mom] = P[ncells-1][.rho]*P[ncells-1][.u]*P[ncells-1][.u] + P[ncells-1][.p]
	F[ncells][.energy] = P[ncells-1][.u]*(P[ncells-1][.rho]*(P[ncells-1][.e] + 0.5*P[ncells-1][.u]*P[ncells-1][.u]) + P[ncells-1][.p])


	for n in 0..<nsteps {
		fmt.printfln("step= %d, t= %.3e", n, f64(n+1)*dt)
		for i in 1..<ncells {
			for cq in Conserved_Quantities {
				F[i][cq] = 0.0
			}
		}
		// Compute all fluxes
		for i in 0..<ncells {
			f_p[i], f_m[i] = vl_flux(P[i])
		}
		// Compute all slopes
		for i in 1..<ncells-1 {
			for cq in Conserved_Quantities {
				s_l := (f_p[i][cq] - f_p[i-1][cq])/dx
				s_r := (f_p[i+1][cq] - f_p[i][cq])/dx
				s_p[i][cq] = flux_slope(s_l, s_r)
				s_l = (f_m[i][cq] - f_m[i-1][cq])/dx
				s_r = (f_m[i+1][cq] - f_m[i][cq])/dx
				s_m[i][cq] = flux_slope(s_l, s_r)
			}
		}
		// Assemble reconstructed flux
		for i in 1..<ncells {
			x := f64(i)*dx
			for cq in Conserved_Quantities {
				F[i][cq] = f_p[i-1][cq] + s_p[i-1][cq]*(x - xs[i-1]) + f_m[i][cq] + s_m[i][cq]*(x - xs[i])
			}
		}
		
		// Compute update
		lambda := dt/dx
		for i in 0..<ncells {
			for cq in Conserved_Quantities {
				U[i][cq] -= lambda*(F[i+1][cq] - F[i][cq])
			}
		}
		// Get primitive quantities up-to-date
		for i in 0..<ncells {
			P[i][.rho] = U[i][.mass]
			P[i][.u] = U[i][.mom]/P[i][.rho]
			P[i][.e] = U[i][.energy]/P[i][.rho] - 0.5*P[i][.u]*P[i][.u]
			P[i][.p] = (GAMMA - 1.0)*P[i][.rho]*P[i][.e]
		}
	}

	f, err := os.open(FILENAME, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o644)
	defer os.close(f)

	fmt.fprintln(f, "x rho p e u")
	for i in 0..<ncells {
		fmt.fprintfln(f, "%.6e %.6e %.6e %.6e %.6e", xs[i], P[i][.rho], P[i][.p], P[i][.e], P[i][.u])
	}
	return
}
