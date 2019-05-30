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
    std::cout << *tree_ << std::endl;
    ComputeForces();
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
        std::cout << "Done interacting with " << *leaf << std::endl
            << "Total force applied: (" << leaf->GetForce(0) << "," << leaf->GetForce(1) << ")" << std::endl;
        // Fetch next leaf and reset the interaction pointer.
        leaf = tree_->GetNextLeaf(leaf);
        node = tree_->GetRoot();
    };
    return E_SUCCESS;
}
