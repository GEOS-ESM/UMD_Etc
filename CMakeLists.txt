option (BUILD_MIT_OCEAN "BUILD_MIT_OCEAN is turned off by default" OFF) #OFF by default
if ( NOT BUILD_MIT_OCEAN )
esma_add_subdirectories(
   UMD_oletkf
   UMD_utils
  )
else()
esma_add_subdirectories(
#   UMD_oletkf
   UMD_utils
  )
endif()

# NOTE There is no real "make" to do in UMD_rc, scripts, or
#      util_hindcast. They are sets of scripts and files. 
#      For now we just INSTALL

# Using MESSAGE_NEVER to not overwhelm the installation
# messages with all these installs
install (
   DIRECTORY UMD_rc
   DESTINATION UMD/rc
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
install (
   DIRECTORY scripts
   DESTINATION UMD/scripts
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
install (
   DIRECTORY scripts
   DESTINATION UMD/scripts
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
install (
   DIRECTORY util_hindcast
   DESTINATION UMD/util_hindcast
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
