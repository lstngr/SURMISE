/** \file stypes.hpp
 * @brief Defines the variables, structures and classes used by the simulation
 * routine.
 */
#ifndef SURMISE_STYPES_
#define SURMISE_STPYES_

#include <string>

/** @brief Possible settings for distributions types (may apply to particles
 * locations, velocities or masses).
 */
enum SDistributions {
    uniform,
    gauss,
    poisson
};

/** @brief Structure containing the simulation parameters.
 * @details Will be read by the root MPI process and passed to other instances
 * before running the simulation.
 */
struct SConfig {
    double dsize;
    unsigned int npart;
    double dt;
    double epsilon;
    SDistributions distr_xlocs;
    SDistributions distr_xvels;
    SDistributions distr_ylocs;
    SDistributions distr_yvels;
    SDistributions distr_mass;
    double unifo_xmin;
    double unifo_xmax;
    double unifo_ymin;
    double unifo_ymax;
    double gauss_xmean;
    double gauss_xstdd;
    double gauss_ymean;
    double gauss_ystdd;
    double unifo_vxmin;
    double unifo_vxmax;
    double unifo_vymin;
    double unifo_vymax;
    double gauss_vxmean;
    double gauss_vxstdd;
    double gauss_vymean;
    double gauss_vystdd;
    double unifo_mass_min;
    double unifo_mass_max;
    double poiss_mass_mean;
    unsigned int max_iter;
    double max_wtime;
    double extra_time;
};

#endif // SURMISE_STYPES
