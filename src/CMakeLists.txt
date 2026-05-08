cmake_minimum_required(VERSION 3.10)
project(Mora LANGUAGES Fortran)

add_executable(mora
    types.f90
    share.f90
    solver.f90
    param.f90
    main.f90
)

set_target_properties(mora PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules"
)

target_include_directories(mora PUBLIC
    "${PROJECT_BINARY_DIR}/modules"
)
