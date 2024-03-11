set(GMP_PREFIX "" CACHE PATH "path ")

find_path(GMP_INCLUDE_DIR gmp.h gmpxx.h 
    PATHS ${GMP_PREFIX}/include /usr/include /usr/local/include )

find_library(GMP_LIBRARY NAMES gmp libgmp 
    PATHS ${GMP_PREFIX}/lib /usr/lib /usr/local/lib)


if(GMP_INCLUDE_DIR AND GMP_LIBRARY)
    get_filename_component(GMP_LIBRARY_DIR ${GMP_LIBRARY} PATH)
    set(GMP_FOUND TRUE)
endif()

if(GMP_FOUND)
   if(NOT GMP_FIND_QUIETLY)
      MESSAGE(STATUS "Found GMP: ${GMP_LIBRARY}")
   endif()
elseif(GMP_FOUND)
   if(GMP_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find GMP")
   endif()
endif()

set(GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})
set(GMP_LIBRARIES ${GMP_LIBRARY})

if(NOT GMP_FOUND)
    # Unset the variables if GMP is not found
    unset(GMP_INCLUDE_DIRS CACHE)
    unset(GMP_LIBRARIES CACHE)
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARIES)