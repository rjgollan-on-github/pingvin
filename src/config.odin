package pingvin

import "core:strings"

import lua "vendor:lua/5.4"

Config :: struct {
    grid2d_file : string,
}

defaults :: `
grid2d_file = ""
`


read_config_from_lua_file :: proc (filename: string) -> (cfg: Config) {
    L := lua.L_newstate()
    lua.open_base(L)
    lua.L_dostring(L, defaults)
    lua.L_dofile(L, strings.unsafe_string_to_cstring(filename))

    lua.getglobal(L, "grid2d_file")
    cfg.grid2d_file = string(lua.tostring(L, -1))
    lua.pop(L, 1)

    return cfg
}


