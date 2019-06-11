#include <stddef.h>
#include <cmath>
#include "stypes.hpp"
#include "serrors.hpp"

template<typename T>
void destroy_vector(std::vector<T*> &v)
{
    while(!v.empty()) {
        delete v.back();
        v.pop_back();
    }
}

MPI_Datatype MPI_Particle;

/** @brief Declares custom MPI datatypes.
 * @details The MPI library provides custom formats for standard C++ types that
 * can be provided to the communication functions (e.g. MPI_Send, MPI_Recv) to
 * indicate how to parse the contents of the provided buffers. In order to use
 * the MPI library's features efficiently, we declare MPI_Datatype objects that
 * allow us to send custom structures, declared within the SURMISE code,
 * efficiently.
 * This function declares a MPI_Particle object to the library by commiting it.
 * If the size of the initially declared MPI_Datatype and the one of the object
 * stored in memory differ, MPI_Type_resize_struct is called.
 */
void make_mpi_types() {
    /* ---------------------------- */
    // Define the MPI Particle type //
    /* ---------------------------- */
    // Inspired by the example 4.17 of the MPI3.1 manual
    MPI_Datatype MPI_Particle_NR;
    // Define structure types and block lengths
    MPI_Datatype ptypes[5]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_LONG};
    int blocklen[5]={2,2,2,1,1};
    // Compute displacements
    MPI_Aint disp[5];
    disp[0] = offsetof( Particle, pos );
    disp[1] = offsetof( Particle, vel );
    disp[2] = offsetof( Particle, frc );
    disp[3] = offsetof( Particle, mass);
    disp[4] = offsetof( Particle, id  );
    MPI_Type_create_struct( 5, blocklen, disp, ptypes, &MPI_Particle_NR );
    // NOTE: Doubles and longs have different alignement (double have 8 bytes,
    // while long only four), such that the compiler might add a four byte
    // padding at the end of the struct. We resize it (will also allow to send
    // multiple particles at once).
    int sz;
    // MPI_Aint lb,ext;
    MPI_Type_size( MPI_Particle_NR, &sz );
    // MPI_Type_get_extent( MPI_Particle_NR, &lb, &ext );
    if(sz != sizeof( Particle )) {
        MPI_Type_create_resized( MPI_Particle_NR,
                0,
                (MPI_Aint)sizeof(Particle),
                &MPI_Particle);
    } else {
        MPI_Particle = MPI_Particle_NR;
    }
    MPI_Type_commit( &MPI_Particle );
}

/** This function frees the types declared to the MPI library by the commits
 * performed in make_mpi_types(). It is called before the main returns.
 */
void free_mpi_types() {
    MPI_Type_free( &MPI_Particle );
}

SConfig::~SConfig() {
    destroy_vector(this->parts);
}

/** Copy operator for a particle.
 * It preserves the index of the copied object.
 * @param[in] other Constant reference to the copied object.
 */
Particle::Particle( const Particle& other )
    :Particle()
{
    *this = other;
}

/** Assignement operator for particles.
 * All properties are conserved, except for the index which is set to -1.
 * @param[in] rhs Right-hand side member of the assignement operator.
 * @returns A modified particle.
 */
Particle& Particle::operator=( const Particle& rhs ) {
    this->pos = rhs.pos;
    this->vel = rhs.vel;
    this->frc = rhs.frc;
    this->mass= rhs.mass;
    // The mass we want to simulate should remain unique in memory. Without
    // explicit instructions, the index is set to -1 to flag that this mass is
    // fictious (i.e. may be assigned to a non-leaf node inside the tree)
    this->id  = -1;
    return *this;
}

/** @brief Adds a particle to the one on the left-hand side of the operator.
 * @details The addition is performed as a center of mass computation, meaning
 * that the LHS of the operator will be the center of mass of the system formed
 * by the two masses. This function is helpful when building the tree and
 * updating it.
 * @param[in] other Particle to add to the left-hand side one.
 * @returns A reference to the modified particle.
 */
Particle& Particle::operator+=( const Particle& other ) {
    double total_mass( this->mass + other.mass );
    this->pos[0] = ( this->mass * this->pos[0] + other.mass * other.pos[0] ) / total_mass;
    this->pos[1] = ( this->mass * this->pos[1] + other.mass * other.pos[1] ) / total_mass;
    this->vel[0] = ( this->mass * this->vel[0] + other.mass * other.vel[0] ) / total_mass;
    this->vel[1] = ( this->mass * this->vel[1] + other.mass * other.vel[1] ) / total_mass;
    this->frc[0] += other.frc[0];
    this->frc[1] += other.frc[1];
    this->mass = total_mass;
    return *this;
}

/** Removes a particle from a center of mass.
 * Has the opposite effect of the overloaded operator+=.
 * @param[in] other Particle to remove from the center of mass.
 * @returns A modified particle.
 */
Particle& Particle::operator-=( const Particle& other ) {
    double total_mass( this->mass - other.mass );
    this->pos[0] = ( this->mass * this->pos[0] - other.mass * other.pos[0] ) / total_mass;
    this->pos[1] = ( this->mass * this->pos[1] - other.mass * other.pos[1] ) / total_mass;
    this->vel[0] = ( this->mass * this->vel[0] - other.mass * other.vel[0] ) / total_mass;
    this->vel[1] = ( this->mass * this->vel[1] - other.mass * other.vel[1] ) / total_mass;
    this->frc[0] -= other.frc[0];
    this->frc[1] -= other.frc[1];
    this->mass = total_mass;
    return *this;
}

/** Overloaded printing operator. It is meant for debugging purposes, one should
 * look into IOManager::WriteOutput to obtain coherently structured data.
 */
std::ostream& operator<<( std::ostream& os, const Particle& p ) {
    os << "[PART " << &p << "] ID=" << p.id << "/MASS=" << p.mass << " POS=(" << p.pos[0] << "," << p.pos[1] << ")";
    return os;
}

/** Computes the gravitational forces between two particles and returns them in
 * an array. Note no force is applied on the passed particles, and that the
 * returned force is the one felt by the first argument due to the second
 * argument.
 */
std::array<double,2> pp_force( const Particle& p1, const Particle& p2 ) {
    // NOTE - We compute the square of the distance instead of working with the
    // Node's functions.
    double strength( - p1.mass * p2.mass / ( ( p1.pos[0]-p2.pos[0] )*( p1.pos[0]-p2.pos[0] ) + ( p1.pos[1]-p2.pos[1] )*( p1.pos[1]-p2.pos[1] ) ) );
    std::array<double,2> force{
        strength * (p1.pos[0] - p2.pos[0]),
        strength * (p1.pos[1] - p2.pos[1])
    };
    return force;
}
