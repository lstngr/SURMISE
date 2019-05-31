#include <iostream>
#include "sroutines.hpp"

Simulation::Simulation( IOManager& io )
    :io_(io), conf_( io.GetConfig() )
{}

Simulation::~Simulation() {
    if(tree_!=NULL) {
        delete tree_;
        tree_ = NULL;
    }
}

SError Simulation::Run() {
    std::cout << "SURMISE run begins." << std::endl;
    BuildTree();
    io_.WriteOutput( *this );
    ComputeForces();
    for( unsigned int iter(0); iter<conf_.max_iter; iter++ ) {
        ComputeForces();
    //    TimeEvolution();
    //    UpdateTree();
        io_.WriteOutput( *this );
    }
    return E_SUCCESS;
}

SError Simulation::BuildTree() {
    std::cout << "Starting tree decomposition." << std::endl;
    tree_ = new QuadTree( conf_ );
    tree_->AddParticle( conf_.parts );
    return E_SUCCESS;
}

SError Simulation::ComputeForces() const {
    Node *leaf(tree_->GetRoot()), *node(leaf);
    leaf = tree_->GetNextLeaf(leaf);
    while( leaf != NULL ){
        leaf->ResetForces();
        while( node != NULL ){
            if( leaf==node ) {
                // A node can't interact with itself, skip if encountered
                node = tree_->GetNext(node);
            } else {
                if( node->GetWidth() / distance(*leaf,*node) >= this->conf_.theta
                        and not node->IsLeaf() ) {
                    // The refinement level is too coarse, move the node pointer
                    // down and retry with children.
                    node = tree_->GetDown(node);
                } else {
                    // The refinement level is good enough. Compute force acting
                    // on the leaf node, and move to the next node horizontally.
                    // NOTE - If everything is correctly implemented, the distance
                    // between node and leaf should be zero, and this statement
                    // shouldn't be reached.
                    leaf->Interact( *node );
                    node = tree_->GetNext(node);
                }
            }
        }
        // Fetch next leaf and reset the interaction pointer.
        leaf = tree_->GetNextLeaf(leaf);
        node = tree_->GetRoot();
    };
    return E_SUCCESS;
}

SError Simulation::TimeEvolution() const {
    Node* leaf = tree_->GetNextLeaf(tree_->GetRoot());
    Particle* p( leaf->GetParticle() );
    while( leaf != NULL ) {
        p->pos[0] += p->vel[0] * conf_.dt;
        p->pos[1] += p->vel[1] * conf_.dt;
        p->vel[0] += p->frc[0] / p->mass * conf_.dt;
        p->vel[1] += p->frc[1] / p->mass * conf_.dt;
        leaf = tree_->GetNextLeaf( leaf );
    }
    return E_SUCCESS;
}
