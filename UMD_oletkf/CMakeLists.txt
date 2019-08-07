esma_set_this ()


set (SRCS
   common.f90 common_mpi.f90 common_mtx.f90 mod_sfmt.f90 
   params_letkf.f90 params_model.f90 params_obs.f90 vars_model.f90
   vars_obs.f90 common_letkf.f90 common_mom4.f90 common_mpi_mom4.f90
   common_obs_mom4.f90 letkf_obs.f90 letkf_local.f90
   letkf_tools.f90 gsw_pot_to_insitu.f90 netlib2.f 
  )

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set( CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${BIG_ENDIAN} ${BYTERECLEN} ${EXTENDED_SOURCE} ${ALIGNCOM}")
endif ()

# MAT There is a bug in common_mom4.f90 that cannot handle bounds-checking
#     For now build with release flags even in debug
set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_RELEASE}")

esma_add_library (${this}
  SRCS ${SRCS}
  INCLUDES ${INC_NETCDF}
  DEPENDENCIES ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES}
  )

target_compile_definitions (${this} PUBLIC _LAPACK_ gmao_intf HAVE_ESMF)

ecbuild_add_executable (
   TARGET oceanda.x 
   SOURCES letkf.f90 
   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
   )

ecbuild_add_executable (
   TARGET oceanobsop.x 
   SOURCES obsop.f90 
   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
   )

ecbuild_add_executable (
   TARGET oceanobs_nc2bin.x 
   SOURCES oceanobs_nc2bin.f90 
   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
   )
