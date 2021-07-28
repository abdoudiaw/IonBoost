########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(IB_src  types.f90
            share.f90
            solver.f90
            param.f90
            IonBoost28.f90
)


# Define the executable in terms of the source files
ADD_EXECUTABLE(${IBEXE} ${IB_src})

SET_PROPERTY (TARGET ${IBEXE} PROPERTY LINKER_LANGUAGE Fortran)

