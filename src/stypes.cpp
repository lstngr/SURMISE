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

Particle::Particle( const Particle& other )
    :Particle()
{
    *this = other;
}

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

Particle& Particle::operator+=( const Particle& other ) {
    double total_mass( this->mass + other.mass );
    this->pos[0] = ( this->mass * this->pos[0] + other.mass * this->pos[0] ) / total_mass;
    this->pos[1] = ( this->mass * this->pos[1] + other.mass * this->pos[1] ) / total_mass;
    this->vel[0] = ( this->mass * this->vel[0] + other.mass * this->vel[0] ) / total_mass;
    this->vel[1] = ( this->mass * this->vel[1] + other.mass * this->vel[1] ) / total_mass;
    this->frc[0] += other.frc[0];
    this->frc[1] += other.frc[1];
    this->mass = total_mass;
    return *this;
}

std::ostream& operator<<( std::ostream& os, const Particle& p ) {
    os << "[PART " << &p << "] ID=" << p.id << "/MASS=" << p.mass << " POS=(" << p.pos[0] << "," << p.pos[1] << ")";
    return os;
}

