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
  INCLUDES ${INC_NETCDF}
  DEPENDENCIES ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES}
  )

target_compile_definitions (${this} PUBLIC _LAPACK_ gmao_intf HAVE_ESMF)

foreach (file ocean_recenter anaice2rst ocean_moments ocean_iau rst2socaice)
   ecbuild_add_executable (
      TARGET ${file}.x 
      SOURCES ${file}.f90 
      LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
      )
endforeach ()

if (F2PY_FOUND)
   add_f2py_module(read_merra2_bcs
      SOURCES read_bin.f90
      DESTINATION bin
      INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
      )
   add_dependencies(read_merra2_bcs ${this})
endif (F2PY_FOUND)

file (GLOB python_files *.py)
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
