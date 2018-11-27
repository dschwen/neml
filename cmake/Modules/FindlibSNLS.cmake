# Variables to set
# libSNLS_FOUND -- managed to find the library
# libSNLS_INCLUDE_DIRS -- path to header files
# libSNLS_LIBRARIES -- path to library

find_path(libSNLS_INCLUDE_DIR
      NAMES SNLS_port.h
      PATHS ${SNLS_PATH}/include
)

find_library(libSNLS_LIBRARY
      NAMES snls
      PATHS ${SNLS_PATH}/lib
)

set(libSNLS_FOUND TRUE)
set(libSNLS_INCLUDE_DIRS ${libSNLS_INCLUDE_DIR})
set(libSNLS_LIBRARIES ${libSNLS_LIBRARY})

mark_as_advanced(libSNLS_INCLUDE_DIR libSNLS_LIBRARY)

