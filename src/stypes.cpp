#include <cmath>
#include "stypes.hpp"

template<typename T>
void destroy_vector(std::vector<T*> &v)
{
    while(!v.empty()) {
        delete v.back();
        v.pop_back();
    }
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

std::ostream& operator<<( std::ostream& os, const Particle& p ) {
    os << "[PART " << &p << "] ID=" << p.id << "/MASS=" << p.mass << " POS=(" << p.pos[0] << "," << p.pos[1] << ")";
    return os;
}

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
