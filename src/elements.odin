package pingvin

// The id of a vertex *is* its index into the global array of vertices.
// We give this a distinct integer type so that we can't use it accidentally
// elsewhere to index into other structures.
VtxId :: distinct int

// A Cartesian vector in 3D space.
Vector3 :: distinct [3]f64

// A 2D face defined by end points.
// We use indices into the global vertices array to 
// define a face.
Face2 :: [2]VtxId

// A quad element is defined by four vertices (hopefully, co-planar)
Quad :: [4]VtxId

// A hexahedral element, defined by 8 vertices
Hex :: [8]VtxId

// Element types as defined by the VTK library. 
// Only a subset is used for our quasi-structured grids.
VTKElement :: enum{line=3, quad=9}

