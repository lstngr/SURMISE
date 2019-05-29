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
