// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cstdlib>

#include "Parallel/Printf/Printf.hpp"
#include "Utilities/Blas.hpp"

namespace {
void print_env_variable(const char* env_var) {
  const char* value = std::getenv(env_var);
  Parallel::printf("%s : %s\n", env_var, value == nullptr ? "Unset" : value);
}
}  // namespace

extern "C" {
#ifdef DISABLE_OPENBLAS_MULTITHREADING
// Declaring this ourselves instead of including cblas.h because our `Blas`
// library does not provide include directories, so cblas.h might not be
// available.
void openblas_set_num_threads(int num_threads);
int openblas_get_num_threads();
#endif  // DISABLE_OPENBLAS_MULTITHREADING
}  // extern "C"

void disable_openblas_multithreading() {
#ifdef DISABLE_OPENBLAS_MULTITHREADING
  print_env_variable("OPENBLAS_NUM_THREADS");
  print_env_variable("OMP_NUM_THREADS");
  Parallel::printf("OpenBLAS current threads: %d\n",
                   openblas_get_num_threads());
  Parallel::printf("Disabling openblas multithreading.\n");
  openblas_set_num_threads(1);
  Parallel::printf("OpenBLAS current threads: %d\n",
                   openblas_get_num_threads());
#endif  // DISABLE_OPENBLAS_MULTITHREADING
}
