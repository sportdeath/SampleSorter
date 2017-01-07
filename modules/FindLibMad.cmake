# Defines
# LIBMAD_FOUND - if found
# LIBMAD_INCLUDE_DIR - the include directory
# LIBMAD_LIBRARY - libmad library path

include(FindPackageHandleStandardArgs)

find_path(LIBMAD_INCLUDE_DIR mad.h)
find_library(LIBMAD_LIBRARY mad)

find_package_handle_standard_args(
  LibMad
  DEFAULT_MSG
  LIBMAD_LIBRARY
  LIBMAD_INCLUDE_DIR
)



