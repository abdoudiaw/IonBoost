cmake_minimum_required(VERSION 3.10)
project(IonBoost LANGUAGES Fortran)

add_executable(ionboost
    types.f90
    config.f90
    linear_solver.f90
    electrons.f90
    grid.f90
    ions.f90
    fields.f90
    pusher.f90
    diagnostics.f90
    main.f90
)

set_target_properties(ionboost PROPERTIES
    Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules"
)

target_include_directories(ionboost PUBLIC
    "${PROJECT_BINARY_DIR}/modules"
)
