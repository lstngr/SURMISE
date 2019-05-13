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
