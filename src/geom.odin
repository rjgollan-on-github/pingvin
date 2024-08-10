package pingvin


/*
 * Compute the volume of a hexahedron.
 *
 * This uses Algoritm (14) in Grandy (1997).
 *
 * *NOTE* Grady labels his vertices different from
 * the choice made here (which follows VTK). As such,
 * the following indices are swapped in his formula:
 *
 * 2 <-> 3
 * 6 <-> 7
 *
 * This implements:
 * 6v_{LD} = [(x6 - x0), (x1 - x0), (x2 - x5)] +
 *           [(x6 - x0), (x4 - x0), (x5 - x7)] +
 *           [(x6 - x0), (x3 - x0), (x7 - x2)]
 *
 * Referece:
 * Grandy (1997)
 * Efficient Computation of Volume of Hexahedral Cells
 * Report no. UCRL-ID-128886, Lawrence Livermore National Laboratory
 *
 * Inputs:
 * hex : a hexahedron defined by its vertices
 *
 * Returns:
 * volume: the volume of the hexahedron
 */
hex_volume :: proc (h: Hex) -> (volume: f64) {
	// A = x6 - x0
	A := global_data.vertices[h[6]] - global_data.vertices[h[0]]
	// B0 = x1 - x0
	B0 := global_data.vertices[h[1]] - global_data.vertices[h[0]]
	// C0 = x2 - x5
	C0 := global_data.vertices[h[2]] - global_data.vertices[h[5]]
	// B1 = x4 - x0
	B1 := global_data.vertices[h[4]] - global_data.vertices[h[0]]
	// C1 = x5 - x7
	C1 := global_data.vertices[h[5]] - global_data.vertices[h[7]]
	// B2 = x3 - x0
	B2 := global_data.vertices[h[3]] - global_data.vertices[h[0]]
	// C2 = x7 - x2
	C2 := global_data.vertices[h[7]] - global_data.vertices[h[2]]

	sum0 := B0.y*C0.z - B0.z*C0.y
	sum0 = sum0 + B1.y*C1.z - B1.z*C1.y
	sum0 = sum0 + B2.y*C2.z - B2.z*C2.y

	sum1 := B0.x*C0.z - B0.z*C0.x
	sum1 = sum1 + B1.x*C1.z - B1.z*C1.x
	sum1 = sum1 + B2.x*C2.z - B2.z*C2.x
	
	sum2 := B0.x*C0.y - B0.y*C0.x
	sum2 = sum2 + B1.x*C1.y - B1.y*C1.x
	sum2 = sum2 + B2.x*C2.y - B2.y*C2.x

	volume = A.x*sum0 - A.y*sum1 + A.z*sum2
	volume = volume/6.0
	return volume
}

bezier_eval :: proc (bez: CubicBezier, t: f64) -> Vector3 {
  /* use de Casteljau algorithm, but unrolled */
  b0 := (1-t)*bez.b0 + t*bez.b1
  b1 := (1-t)*bez.b1 + t*bez.b2
  b2 := (1-t)*bez.b2 + t*bez.b3
  //
  b0 = (1-t)*b0 + t*b1
  b1 = (1-t)*b1 + t*b2
  //
  b0 = (1-t)*b0 + t*b1
  return b0
}

bezier_t_from_x :: proc (bez: CubicBezier, x: f64) -> (t: f64) {
	t = (x - bez.b0.x)/(bez.b3.x - bez.b0.x)
	return
}


