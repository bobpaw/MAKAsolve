find_path(AmgX_INCLUDE_DIR amgx_c.h
  HINTS ${AmgX_DIR} $ENV{AmgX_DIR}
  PATH_SUFFIXES include
)

find_library(AmgX_LIBRARY
  NAMES amgx
  HINTS ${AmgX_DIR} $ENV{AmgX_DIR}
  PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AmgX DEFAULT_MSG
  AmgX_LIBRARY AmgX_INCLUDE_DIR
)

if(AmgX_FOUND)
  add_library(amgx UNKNOWN IMPORTED)
  set_target_properties(amgx PROPERTIES
    IMPORTED_LOCATION "${AmgX_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${AmgX_INCLUDE_DIR}"
  )
endif()
