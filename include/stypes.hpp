/** \file stypes.hpp
 * @brief Defines the variables, structures and classes used by the simulation
 * routine.
 */
#ifndef SURMISE_STYPES_
#define SURMISE_STYPES_

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
    /** Configuration destructor. Assumes masses were initialized dynamically be
     * the "usual" routine (and not the user's nasty practices).
     * @todo Make the call safer, or protect the particle storage.
     */
    ~SConfig(){ delete[] parts; }
    double dsize; //!< Domain size (square region's edge)
    unsigned int npart; //!< Number of Particles
    double dt; //!< Time step for temporal evolution.
    double epsilon; //!< Precision for the Fast Multipole Method
    unsigned int max_iter; //!< Maximum number of iterations.
    double max_wtime; //!< Maximum wall time allowed.
    double extra_time; //!< Extra time (substracted from max_wtime), allows exit routines to be performed peacefully.
    Particle *parts; //!< Particles to be simulated.
};

#endif // SURMISE_STYPES
