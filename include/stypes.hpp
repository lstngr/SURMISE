/** \file stypes.hpp
 * @brief Defines the variables, structures and classes used by the simulation
 * routine.
 */
#ifndef SURMISE_STYPES_
#define SURMISE_STPYES_

#include <string>

/** Describes a particle in space. Contains relevant informations to
 * itself (position, velocity, mass) as well as flags set by the domain classes.
 */
struct Particle {
    float pos[2]; //!< 2D Position of the particle
    float vel[2]; //!< 2D Velocity of the particle
    float mass; //!< Mass of the particle
    int id; //!< Mass identifier (so we do not mix them during time-evolution).
};

/** @brief Structure containing the simulation parameters.
 * @details Will be read by the root MPI process and passed to other instances
 * before running the simulation.
 */
struct SConfig {
    ~SConfig(){ delete[] parts; }
    double dsize;
    unsigned int npart;
    double dt;
    double epsilon;
    unsigned int max_iter;
    double max_wtime;
    double extra_time;
    Particle *parts;
};

#endif // SURMISE_STYPES
