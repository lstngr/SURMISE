/** \file stypes.hpp
 * @brief Defines the variables, structures and classes used by the simulation
 * routine.
 */
#ifndef SURMISE_STYPES_
#define SURMISE_STYPES_

#include <string>
#include <vector>
#include <array>
#include <iostream>

/** Describes a particle in space. Contains relevant informations to
 * itself (position, velocity, mass) as well as flags set by the domain classes.
 */
struct Particle {
    std::array<double,2> pos; //!< 2D Position of the particle
    std::array<double,2> vel; //!< 2D Velocity of the particle
    std::array<double,2> frc; //!< 2D Force acting on the particle
    double mass; //!< Mass of the particle
    long id; //!< Mass identifier (so we do not mix them during time-evolution).
    Particle(){}
    Particle( const Particle& other );
    Particle& operator=( const Particle& rhs );
    Particle& operator+=( const Particle& other );
};
Particle operator+( Particle lhs, const Particle& rhs );
std::ostream& operator<<( std::ostream& os, const Particle& p );

/** @brief Structure containing the simulation parameters.
 * @details Will be read by the root MPI process and passed to other instances
 * before running the simulation.
 */
struct SConfig {
    /** Configuration destructor. Assumes masses were initialized dynamically be
     * the "usual" routine (and not the user's nasty practices).
     * @todo Make the call safer, or protect the particle storage.
     */
    ~SConfig();
    double dsize; //!< Domain size (square region's edge)
    unsigned int npart; //!< Number of Particles
    double dt; //!< Time step for temporal evolution.
    double epsilon; //!< Precision for the Fast Multipole Method
    unsigned int max_iter; //!< Maximum number of iterations.
    double max_wtime; //!< Maximum wall time allowed.
    double extra_time; //!< Extra time (substracted from max_wtime), allows exit routines to be performed peacefully.
    std::vector<Particle*> parts; //!< Particles to be simulated.
};

#endif // SURMISE_STYPES
