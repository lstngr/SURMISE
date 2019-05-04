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
    double yo; //!< Amount of coolness of the simulation
    std::string str; //!< Punchline of the simulator
};

#endif // SURMISE_STYPES
