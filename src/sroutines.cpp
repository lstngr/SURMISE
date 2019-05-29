#include <iostream>
#include "sroutines.hpp"

Simulation::~Simulation() {
    if(tree_!=NULL)
        delete tree_;
}

SError Simulation::Run() {
    std::cout << "Coucou" << std::endl;
    return E_SUCCESS;
}
