#include <iostream>
#include "sroutines.hpp"

Simulation::~Simulation() {
    if(tree_!=NULL) {
        delete tree_;
        tree_ = NULL;
    }
}

SError Simulation::Run() {
    std::cout << "SURMISE run begins." << std::endl;
    BuildTree();
    // for( unsigned int iter(0); iter<conf_.max_iter; iter++ ) {
    //     ComputeForces();
    //     TimeEvolution();
    //     UpdateTree();
    // }
    return E_SUCCESS;
}

SError Simulation::BuildTree() {
    std::cout << "Starting tree decomposition." << std::endl;
    tree_ = new QuadTree( conf_ );
    tree_->AddParticle( conf_.parts );
    return E_SUCCESS;
}
