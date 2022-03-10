esma_set_this ()

esma_add_subdirectories(sosie-2.6.4)

set (SRCS
   ocean_recenter_utils.f90 ocnice_utils.f90
   )

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set( CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${BIG_ENDIAN} ${BYTERECLEN} ${EXTENDED_SOURCE} ${ALIGNCOM}")
endif ()

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES NetCDF::NetCDF_Fortran MPI::MPI_Fortran
  )

target_compile_definitions (${this} PUBLIC _LAPACK_ gmao_intf HAVE_ESMF)

foreach (file ocean_recenter anaice2rst ocean_moments ocean_iau)
   ecbuild_add_executable (
      TARGET ${file}.x
      SOURCES ${file}.f90
      LIBS ${this} NetCDF::NetCDF_Fortran MPI::MPI_Fortran Threads::Threads
      )
endforeach ()

if (USE_F2PY)
   find_package(F2PY2)
   if (F2PY2_FOUND)
      esma_add_f2py2_module(read_merra2_bcs
         SOURCES read_bin.f90
         DESTINATION bin
         INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
         )
      add_dependencies(read_merra2_bcs ${this})
   endif ()
endif ()

file (GLOB python_files CONFIGURE_DEPENDS *.py)
install(
   PROGRAMS ${python_files}
   DESTINATION bin
   )

install (
   DIRECTORY plot
   DESTINATION UMD/utils/plot
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
