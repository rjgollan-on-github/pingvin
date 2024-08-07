package pingvin

import "core:testing"
import "core:fmt"

@(test)
hex_volume_test :: proc (t: ^testing.T) {
	global_data.vertices = make([dynamic]Vector3, 8)
	defer delete(global_data.vertices)
	v := global_data.vertices[:]
	v[0] = {0.0, 0.0, 0.0}
	v[1] = {1.0, 0.0, 0.0}
	v[2] = {1.0, 1.0, 0.0}
	v[3] = {0.0, 1.0, 0.0}
	v[4] = {0.0, 0.0, 1.0}
	v[5] = {1.0, 0.0, 1.0}
	v[6] = {1.0, 1.0, 1.0}
	v[7] = {0.0, 1.0, 1.0}

	hex := Hex{0, 1, 2, 3, 4, 5, 6, 7}
	vol := hex_volume(hex)

	expected_volume := 1.0
	testing.expectf(t, abs(vol - expected_volume) <= 1.0e-6, "hex= %v, vol= %v", hex, vol)

	// Now, reduce z values by half
	for &vv in v {
		vv.z = 0.5*vv.z
	}
	vol = hex_volume(hex)
	expected_volume = 0.5
	testing.expectf(t, abs(vol - expected_volume) <= 1.0e-6, "hex= %v, vol= %v", hex, vol)

	// and reduce x and y by half
	for &vv in v {
		vv.x = 0.5*vv.x
		vv.y = 0.5*vv.y
	}
	vol = hex_volume(hex)
	expected_volume = 0.125
	testing.expectf(t, abs(vol - expected_volume) <= 1.0e-6, "hex= %v, vol= %v", hex, vol)
}
