package pingvin

FileOpenFailed :: struct {
    filename: string,
}

DirectoryOpenFailed :: struct {
    directory: string,
}

GridReadingError :: struct {
    msg: string,
}

GridFormationError :: struct {
    msg: string,
}
