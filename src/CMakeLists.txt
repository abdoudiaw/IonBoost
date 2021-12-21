# CMake project file for Ionboost

cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(Ionboost VERSION 28.0)

ENABLE_LANGUAGE( Fortran )


add_executable(lmp types.f90
                        share.f90
                        solver.f90
                        param.f90
                        IonBoost28.f90)

target_include_directories(lmp PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )
