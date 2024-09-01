package pingvin

import "core:strings"
import "core:fmt"
import "core:os"
import "core:log"

import lua "vendor:lua/5.4"

Config :: struct {
    grid2d_file :                string,
    cross_section_dir :          string,
    output_vtk_file :            string,
    n_xsects :                   int,
    dx :                         f64,
    Mach_inflow :                f64,
    p_inflow :                   f64,
    T_inflow :                   f64,
    // Solver settings
    max_newton_steps:            int,
    slice_relative_residual:     f64,
    slice_change_in_update:      f64,
    max_gmres_iterations:        int,
    perturbation_size:           f64,
    gmres_relative_residual:     f64,
}

defaults :: `
grid2d_file = ""
cross_section_dir = "xsect"
output_vtk_file = "pingvin-flow-field.vtk"
dx = -1.0
max_newton_steps = 10
slice_relative_residual = 1.0e-6
slice_change_in_update = 1.0e-6
max_gmres_iterations = 10
perturbation_size = 1.0e-250
gmres_relative_residual = 1.0e-3
`

lua_get_optional_string :: proc (L: ^lua.State, field: cstring) -> (result: string, found: bool) {
    found = false
    result = ""
    lua.getglobal(L, field)
    if !lua.isnil(L, -1) {
        result = string(lua.tostring(L, -1))
        found = true
    }
    lua.pop(L, 1)
    return result, found
}

lua_get_integer :: proc (L: ^lua.State, field : cstring, err_msg: string) -> (result: int) {
    lua.getglobal(L, field)
    if !lua.isnil(L, -1) {
        result = int(lua.tointeger(L, -1))
        lua.pop(L, 1)
        return result
    }
    // If we get here, there was an error finding the field
    fmt.println(err_msg)
    os.exit(1)
}

lua_get_optional_integer :: proc (L: ^lua.State, field: cstring) -> (result: int, found: bool) {
    found = false
    lua.getglobal(L, field)
    if !lua.isnil(L, -1) {
        result = int(lua.tointeger(L, -1))
        found = true
    }
    lua.pop(L, 1)
    return result, found
}

lua_get_float :: proc (L: ^lua.State, field : cstring, err_msg: string) -> (result: f64) {
    lua.getglobal(L, field)
    if !lua.isnil(L, -1) {
        result = f64(lua.tonumber(L, -1))
        lua.pop(L, 1)
        return result
    }
    // If we get here, there was an error finding the field
    fmt.println(err_msg)
    os.exit(1)
}

lua_get_optional_float :: proc (L: ^lua.State, field: cstring) -> (result: f64, found: bool) {
    found = false
    lua.getglobal(L, field)
    if !lua.isnil(L, -1) {
        result = f64(lua.tonumber(L, -1))
        found = true
    }
    lua.pop(L, 1)
    return result, found
}

read_config_from_lua_file :: proc (filename: string) -> (cfg: Config) {
    L := lua.L_newstate()
    lua.open_base(L)
    lua.L_dostring(L, defaults)
    lua.L_dofile(L, strings.unsafe_string_to_cstring(filename))

    found : bool
    str_result : string
    int_result : int
    float_result : f64

    str_result, found = lua_get_optional_string(L, "grid2d_file")
    if found do cfg.grid2d_file = str_result

    str_result, found = lua_get_optional_string(L, "cross_section_dir")
    if found do cfg.cross_section_dir = str_result

    str_result, found = lua_get_optional_string(L, "output_vtk_file")
    if found do cfg.output_vtk_file = str_result

    err_msg := "PVN-ERROR: 'no_cross_sections' not set in job input file. An integer is expected."
    cfg.n_xsects = lua_get_integer(L, "no_cross_sections", err_msg)

    err_msg = "PVN-ERROR: 'dx' not set in job input file. A floating point value is expected."
    cfg.dx = lua_get_float(L, "dx", err_msg)
    if (cfg.dx <= 0.0) {
        fmt.printfln("PVN-ERROR: Error in setting 'dx'. Value must be positive. Value is %v", cfg.dx)
    }

    err_msg = "PVN-ERROR: 'Mach_inflow' not set in job input file. A floating point value is expected."
    cfg.Mach_inflow = lua_get_float(L, "Mach_inflow", err_msg)

    err_msg = "PVN-ERROR: 'p_inflow' not set in job input file. A floating point value is expected."
    cfg.p_inflow = lua_get_float(L, "p_inflow", err_msg)
    

    err_msg = "PVN-ERROR: 'T_inflow' not set in job input file. A floating point value is expected."
    cfg.T_inflow = lua_get_float(L, "T_inflow", err_msg)
    
    // Parameters related to solver settings
    int_result, found = lua_get_optional_integer(L, "max_newton_steps")
    if found do cfg.max_newton_steps = int_result

    float_result, found = lua_get_optional_float(L, "slice_relative_residual")
    if found do cfg.slice_relative_residual = float_result

    float_result, found = lua_get_optional_float(L, "slice_change_in_update")
    if found do cfg.slice_change_in_update = float_result

    int_result, found = lua_get_optional_integer(L, "max_gmres_iterations")
    if found do cfg.max_gmres_iterations = int_result

    float_result, found = lua_get_optional_float(L, "perturbation_size")
    if found do cfg.perturbation_size = float_result

    float_result, found = lua_get_optional_float(L, "gmres_relative_residual")
    if found do cfg.gmres_relative_residual = float_result

    return cfg
}


