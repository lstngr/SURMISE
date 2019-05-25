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

std::array<double,2> Particle::PForce( const Particle* p ) const {
    double dist2( (this->pos[0]-p->pos[0])*(this->pos[0]-p->pos[0]) + (this->pos[1]-p->pos[1])*(this->pos[1]-p->pos[1])  );
    return std::array<double,2>{
        - this->mass * p->mass / dist2 * (this->pos[0] - p->pos[0]),
        - this->mass * p->mass / dist2 * (this->pos[1] - p->pos[1]),
    };
}
