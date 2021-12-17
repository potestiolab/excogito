# ensuring math is present
set(CMAKE_REQUIRED_INCLUDES math.h)
set(CMAKE_REQUIRED_LIBRARIES m)
add_library(lib ini.c alignment.c hierarchical_clustering.c observables.c mapping.c geometry.c io.c my_malloc.c sampling.c traj.c optimize.c random_mapping.c measure.c norm.c cosine.c distance.c measure_kl.c optimize_kl.c random_mapping_kl.c)
