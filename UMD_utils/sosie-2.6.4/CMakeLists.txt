esma_set_this ()

set (SRCS
   src/io_ezcdf.f90
   src/mod_akima_1d.f90
   src/mod_bicub_2d.f90
   src/mod_conf.f90
   src/mod_grids.f90
   src/mod_interp.f90
   src/mod_akima_2d.f90
   src/mod_bilin_2d.f90
   src/mod_drown.f90
   src/mod_init.f90
   src/mod_manip.f90
   )

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set( CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${BIG_ENDIAN} ${BYTERECLEN} ${EXTENDED_SOURCE} ${ALIGNCOM}")
endif ()

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES NetCDF::NetCDF_Fortran MPI::MPI_Fortran
  )

target_compile_definitions (${this} PUBLIC _LAPACK_ gmao_intf HAVE_ESMF)

foreach (file sosie corr_vect interp_to_line mask_drown_field)
   ecbuild_add_executable (
      TARGET ${file}.x
      SOURCES src/${file}.f90
      LIBS ${this} NetCDF::NetCDF_Fortran MPI::MPI_Fortran Threads::Threads
      )
endforeach ()

