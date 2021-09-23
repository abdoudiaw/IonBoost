#Just adding this for test sake.
# CMake project file for KG

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2)
PROJECT(IONBOOST)

ENABLE_LANGUAGE( Fortran )

# Set the KG version
SET(VERSION 1.2.3)


##################################################
#  Define options
##################################################



IF (USE_OPENMP)
   SET(CMAKE_Fortran_FLAGS    "${OpenMP_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS}")
   SET(CMAKE_EXE_LINKER_FLAGS "${OpenMP_Fortran_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
   add_definitions( -Dopen_mp_only )
ENDIF ()




############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(KGEXE mix)

ADD_SUBDIRECTORY(src)

