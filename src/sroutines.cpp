#include <iostream>
#include <cmath>
#include "sroutines.hpp"

/** Instantiates a Simulation object. This class requires an input/output
 * manager to be initialized by the user, as it will fetch its configuration
 * from it.
 * @param[in] io Reference to an initialized IOManager.
 */
Simulation::Simulation( IOManager& io )
    :io_(io), timer_(new Timer()), conf_( io.GetConfig() ), tree_(NULL), update_list_(NULL)
{}

/** Destroys the simulation object by freeing up the QuadTree allocated by the
 * class. All other members should be well behaved (the IOManager is set up in
 * the main, and the configuration object depends on it).
 */
Simulation::~Simulation() {
    if(tree_!=NULL) {
        delete tree_;
        tree_ = NULL;
    }
    if(update_list_!=NULL) {
        delete[] update_list_;
        update_list_=NULL;
    }
    if( timer_ != NULL ) {
        delete timer_;
        timer_=NULL;
    }
}

/** @brief Runs a N-Body simulation using the Barnes-Hut algorithm.
 * @details This routine is responsible for managing a full simulation, once the
 * class is properly initialized. The routine first builds a tree from scratch,
 * and performs the following calls every time step:
 *
 * - ComputeForces: Iterates on all leafs of the tree and computes the forces
 *   acting on their center of mass according to the Barnes-Hut algorithm.
 * - TimeEvolution: Performs a temporal evolution of the system using an Euler
 *   scheme.
 * - UpdateTree: Scans for updates needed by the tree. This includes particle
 *   reassignement to other nodes, node creation and destruction according to
 *   their occupancy, as well as upstream update propagation to keep the center
 *   of mass in higher levels up to date.
 * - IOManager::WriteOutput: Calls the output routine as defined in the
 *   IOManager class. It outputs leaf information, tree structure, etc.
 *
 * Once the required number of iterations has been performed, the routine
 * returns an error code indicating the success or failure of the algorithm.
 *
 * @returns An error code.
 */
SError Simulation::Run() {
    int rank,size;
    double start_time(MPI_Wtime())
        ,wtime(0.);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    if( rank==0 ){
        wtime = MPI_Wtime();
        std::cout << "["<< wtime << "s] SURMISE run begins." << std::endl
            << "   MPI_WTIME_IS_GLOBAL=" << (bool)MPI_WTIME_IS_GLOBAL << std::endl;
    }
    timer_->StartTimer( T_DISTR );
    io_.DistributeParticles(*this);
    timer_->StopTimer( T_DISTR );
    for( unsigned iter(0); iter<conf_.max_iter; iter++ ) {
        timer_->StartTimer( T_ITER );
        ComputeForces();
        TimeEvolution();
        // Sync Leafs, includes gathering to master
        timer_->StartTimer( T_LEAFSYNC );
        io_.SyncLeafs(*this);
        timer_->StopTimer( T_LEAFSYNC );
        UpdateTree();
        timer_->StartTimer( T_OUTPUT );
        if(rank==0 and iter%conf_.write_freq ){
            io_.WriteOutput( *this );
        } else if( iter%conf_.write_freq ){
            double *sbuf(new double[T_COUNT]), *rbuf(NULL);
            for( unsigned i_timer(0); i_timer<T_COUNT; i_timer++ )
                sbuf[i_timer] = timer_->GetTime(i_timer);
            MPI_Gather( sbuf, T_COUNT, MPI_DOUBLE, rbuf, T_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD );
            delete[] sbuf;
        }
        timer_->StopTimer( T_OUTPUT );
        timer_->StopTimer( T_ITER );
        // If sufficient unbalancing, recompute indicies.
        // If time elapsed, cleanfully exit
        MPI_Barrier( MPI_COMM_WORLD );
        if( MPI_Wtime()-start_time > conf_.max_wtime - conf_.extra_time ){
            if( rank==0 )
                std::cout << "[" << MPI_Wtime() << "s] Run time elapsed, aborting"
                    << " further computation at iter=" << iter << "." << std::endl;
            return E_TIMEOUT;
        }
    }
    return E_SUCCESS;
}

/** Allocates a tree and sends all particles found in the configuration object
 * to it.
 * @returns An error code.
 */
SError Simulation::BuildTree() {
    tree_ = new QuadTree( conf_ );
    tree_->AddParticle( conf_.parts );
    return E_SUCCESS;
}

/** @brief Clean-up routine for the tree.
 * @details This method cleans the tree at each iteration. It is responsible
 * for:
 *
 * - Moving particles that exited their Node's coverage area to a proper node,
 *   and allocate the subsequent node if required.
 * - Deleting nodes that have been emptied of their particles, and eventually
 *   update the local tree structure if we find unnecessarily deep levels there.
 *   This might happen when a parent node posses two leafs, and one of the
 *   leafs moves out of the parent's coverage area: the remaining leaf can be
 *   moved on level up and replace the parent node. This is replacement can
 *   happen as many time as unnecessary parent nodes are found.
 * - Clearing and updating upstream centers of mass. After a time evolution
 *   step, the upstream nodes informations are out of date compared to the leaf
 *   levels (leafs have moved, and the structure of the tree may have changed).
 *   The upstream levels are first cleared, and the information of each leaf is
 *   propagated upwards in the tree.
 *
 * @returns An error code.
 * @note The upstream center of mass update is not optimized. As the structure
 * of the tree can change after the node pruning step, we cannot perform it as
 * the leafs are updated in a straightforward manner. For this reason, this part
 * of the routine is implemented in a "brute-force" manner.
 * @warning If this method is called when only one particle remains in the
 * simulation, a segmentation fault occurs within this function.
 */
SError Simulation::UpdateTree() const {
    timer_->StartTimer( T_TREEUPDATE );
    Node *leaf( tree_->GetNextLeaf( tree_->GetRoot() ) );
    // Reassignement step. Particles are moved where they're supposed to be.
    do {
        // If the particle moved out of its node's coverage area, reassignement
        // is needed.
        // The Add-/RemoveParticle methods update upstream informations such as
        // center of mass and position.
        Particle *p ( leaf->GetParticle() );
        if( not leaf->BelongsTo(p) ) {
            tree_->RemoveParticle( leaf );
            tree_->AddParticle( p );
            // The node is empty now. We first move the leaf pointer before
            // changing the tree structure, and prune the node.
            leaf = tree_->GetNextLeaf( tree_->PruneNode( leaf ) );
            // Local tree structure may change on several levels. Safer to
            // restart from Root. TODO(Might improve? Fetch parents, go down
            // again?)
            // leaf = tree_->GetNextLeaf( tree_->GetRoot() );
        } else {
            leaf = tree_->GetNextLeaf( leaf );
        }
    } while( leaf != NULL );
    // Upstream center of mass update. All the center of mass are updated
    // according to leaf nodes. A first pass resets them, a second one updates
    // them with leaf info.
    leaf = tree_->GetNextLeaf( tree_->GetRoot() );
    if( leaf==NULL )
        return E_SUCCESS;
    do {
        Node* up( leaf->GetParent() );
        do {
            up->ResetCenterOfMass();
            up = up->GetParent();
        } while( up != NULL );
        leaf = tree_->GetNextLeaf(leaf);
    } while( leaf != NULL );
    // Recursively add leaf level influence.
    leaf = tree_->GetNextLeaf( tree_->GetRoot() );
    do {
        Node* up( leaf->GetParent() );
        do {
            *(up->GetParticle()) += *(leaf->GetParticle());
            up = up->GetParent();
        } while( up != NULL );
        leaf = tree_->GetNextLeaf(leaf);
    } while( leaf != NULL );
    timer_->StopTimer( T_TREEUPDATE );
    return E_SUCCESS;
}

/** Iterates over all leafs in the tree, first clears their force attribute, and
 * computes a new force based on the Barnes-Hut algorithm.
 * @returns An error code.
 * @note The Barnes-Hut criterion for computing an interaction with a distant
 * node reads s/d < theta. We implement the criterion as s^2 < theta^2 * d^2
 * since benchmarks show an performance improvement of 10% if a call to
 * std::sqrt and a double division is avoided.
 */
SError Simulation::ComputeForces() const {
    timer_->StartTimer( T_FORCES );
    Node *leaf(tree_->GetRoot()), *node(leaf);
    leaf = tree_->GetNextLeaf(leaf);
    while( leaf != NULL ){
        leaf->ResetForces();
        if( update_list_[leaf->GetParticle()->id] ) {
            while( node != NULL ){
                if( leaf==node ) {
                    // A node can't interact with itself, skip if encountered
                    node = tree_->GetNext(node);
                } else {
                    if( node->GetWidth()*node->GetWidth() >= this->conf_.theta *
                            this->conf_.theta * distance2(*leaf,*node) and not
                            node->IsLeaf() ) {
                        // The refinement level is too coarse, move the node
                        // pointer down and retry with children.
                        node = tree_->GetDown(node);
                    } else {
                        // The refinement level is good enough. Compute force
                        // acting on the leaf node, and move to the next node
                        // horizontally.
                        // NOTE - If everything is correctly implemented, the
                        // distance between node and leaf should be zero, and
                        // this statement shouldn't be reached.
                        leaf->Interact( *node );
                        node = tree_->GetNext(node);
                    }
                }
            }
        }
        // Fetch next leaf and reset the interaction pointer.
        leaf = tree_->GetNextLeaf(leaf);
        node = tree_->GetRoot();
    }
    timer_->StopTimer( T_FORCES );
    return E_SUCCESS;
}

/** Applies an Euler scheme to all leafs of the tree, assuming the forces acting
 * on the particles have been computed beforehand.
 * @returns An error code.
 */
SError Simulation::TimeEvolution() const {
    timer_->StartTimer( T_EVOL );
    Node* leaf = tree_->GetNextLeaf(tree_->GetRoot());
    while( leaf != NULL ) {
        if( update_list_[leaf->GetParticle()->id] ) {
            Particle* p( leaf->GetParticle() );
            p->pos[0] += p->vel[0] * conf_.dt;
            p->pos[1] += p->vel[1] * conf_.dt;
            p->vel[0] += p->frc[0] / p->mass * conf_.dt;
            p->vel[1] += p->frc[1] / p->mass * conf_.dt;
        }
        leaf = tree_->GetNextLeaf( leaf );
    }
    timer_->StopTimer( T_EVOL );
    return E_SUCCESS;
}
