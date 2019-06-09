#ifndef SURMISE_STYPES_
#define SURMISE_STYPES_

#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <mpi.h>

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
    Particle& operator-=( const Particle& other );
};
std::ostream& operator<<( std::ostream& os, const Particle& p );
std::array<double,2> pp_force( const Particle& p1, const Particle& p2 );

/** MPI Datatype corresponding to a Particle structure.
 * \see make_mpi_types
 */
extern MPI_Datatype MPI_Particle;
void make_mpi_types();
void free_mpi_types();

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
    double theta; //!< Barnes-Hut "precision" parameter, controlling how coarse the interaction between distant nodes should be.
    unsigned int max_iter; //!< Maximum number of iterations.
    double max_wtime; //!< Maximum wall time allowed.
    double extra_time; //!< Extra time (substracted from max_wtime), allows exit routines to be performed peacefully.
    unsigned write_freq; //!< Output frequency. The program will output data every `write_freq` time steps starting from the first iteration.
    std::vector<Particle*> parts; //!< Particles to be simulated.
    std::string ipath={}; //!< Input path to the configuration files.
    std::string opath={}; //!< Output path where data should be written.
};

#endif // SURMISE_STYPES
