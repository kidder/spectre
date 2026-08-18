# Distributed under the MIT License.
# See LICENSE.txt for details.

find_package(BLAS REQUIRED)

set(IS_OPENBLAS FALSE)

if(BLAS_FOUND)
  string(TOLOWER "${BLAS_LIBRARIES}" BLAS_LIBS_LOWER)
  if(BLAS_LIBS_LOWER MATCHES "openblas")
    set(IS_OPENBLAS TRUE)
  else()
    if(TARGET BLAS::BLAS)
      get_target_property(BLAS_LINK_LIBS BLAS::BLAS INTERFACE_LINK_LIBRARIES)
      string(TOLOWER "${BLAS_LINK_LIBS}" BLAS_TARGET_LIBS_LOWER)
      if(BLAS_TARGET_LIBS_LOWER MATCHES "openblas")
        set(IS_OPENBLAS TRUE)
      endif()
    endif()
  endif()
endif()

if(IS_OPENBLAS)
  message(STATUS "Detected BLAS vendor: OpenBLAS")
endif()

message(STATUS "BLAS libs: " ${BLAS_LIBRARIES})
file(APPEND
  "${CMAKE_BINARY_DIR}/BuildInfo.txt"
  "BLAS_LIBRARIES: ${BLAS_LIBRARIES}\n"
  )

set_property(
  GLOBAL APPEND PROPERTY SPECTRE_THIRD_PARTY_LIBS
  BLAS::BLAS
  )

# Check if we have found OpenBLAS and can disable its multithreading, since it
# conflicts with Charm++ parallelism. Details:
# https://github.com/xianyi/OpenBLAS/wiki/Faq#multi-threaded
# We use `execute_process` instead of `try_compile` to avoid potentially slow
# disk IO.
set(
  CHECK_DISABLE_OPENBLAS_MULTITHREADING_SOURCE
  "extern \"C\" { void openblas_set_num_threads(int); }\n\
int main() { openblas_set_num_threads(1); }"
  )
string(REPLACE ";" " " BLAS_LIBRARIES_JOINED_WITH_SPACES "${BLAS_LIBRARIES}")
execute_process(
  COMMAND
  bash -c
  "${CMAKE_CXX_COMPILER} ${BLAS_LIBRARIES_JOINED_WITH_SPACES} -x c++ - <<< $'\
${CHECK_DISABLE_OPENBLAS_MULTITHREADING_SOURCE}' -o /dev/null"
  RESULT_VARIABLE CHECK_DISABLE_OPENBLAS_MULTITHREADING_RESULT
  OUTPUT_VARIABLE CHECK_DISABLE_OPENBLAS_MULTITHREADING_OUTPUT
  ERROR_VARIABLE CHECK_DISABLE_OPENBLAS_MULTITHREADING_ERROR
  )
if(${CHECK_DISABLE_OPENBLAS_MULTITHREADING_RESULT} EQUAL 0)
  set(DISABLE_OPENBLAS_MULTITHREADING ON)
  add_definitions(-DDISABLE_OPENBLAS_MULTITHREADING)
  message(STATUS "Disabled OpenBLAS multithreading")
else()
  if(IS_OPENBLAS)
    message(STATUS
            "BLAS vendor is OpenBLAS but disabling multithreading failed.")
    message(STATUS "Output: ${CHECK_DISABLE_OPENBLAS_MULTITHREADING_OUTPUT}")
    message(STATUS "Error: ${CHECK_DISABLE_OPENBLAS_MULTITHREADING_ERROR}")
    set(DISABLE_OPENBLAS_MULTITHREADING ON)
    add_definitions(-DDISABLE_OPENBLAS_MULTITHREADING)
  else()
    message(STATUS "BLAS vendor is not OpenBLAS. Make sure it doesn't "
      "try to do multithreading that might conflict with Charm++ parallelism.")
  endif()
endif()
