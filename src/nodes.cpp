#include <iostream>
#include <iomanip>
#include <cmath>
#include "nodes.hpp"

//////////////////// NODE ///////////////////////

Node::Node( const SConfig& conf )
    :Node(NULL,NULL)
{
    // Level, index, left and bottom are already initialized correctly.
    top = conf.dsize; right = conf.dsize;
}

Node::Node( Node* parent, Particle* part )
    :parent_(parent), children_{NULL,NULL,NULL,NULL}, level_(0), index_(0), com_(part), left(0.0), right(0.0), top(0.0), bottom(0.0)
{
    // If parent exists, we can infer some properties (such as index, level,
    // boundaries) to the current node.
    if(parent_!=NULL) {
        level_ = parent_->level_++;
        // Find which child we are and compute properties
        for( unsigned int ic(0); ic<4; ic++ ){
            if( parent_->children_[ic] == this ){
                // Shift the node's index by two bits compared to the parent,
                // then, add local index
                index_ = (parent_->index_ << 2) + ic;
                // Check child location relative to parent. Using the Z-curve
                // structure, if the zeroth bit of ic is one (ic=0,2), the node
                // is north. If the next bit of ic is one (ic=2,3), the node is
                // east.
                bool is_north( ic & 1U );
                bool is_east( (ic>>1) & 1U );
                left   = parent_->left   + 0.5*(double)(is_east)  *(parent_->right-parent_->left);
                right  = parent_->right  - 0.5*(double)(!is_east) *(parent_->right-parent_->left);
                bottom = parent_->bottom + 0.5*(double)(is_north) *(parent_->top-parent_->bottom);
                top    = parent_->top    - 0.5*(double)(!is_north)*(parent_->top-parent_->bottom);
                // We found the node, break iteration
                break;
            }
        }
    }
}

/** @brief Creates a generic node.
 * @param[in] parent Pointer to a parent node. If none is provided, the pointer
 * is left null and the node is interpreted to be at the root.
 * @param[in] idx Index of the child within its parent's children subtree.
 */
Node::Node( Node *parent, unsigned int idx )
    :parent_(parent), children_{NULL,NULL,NULL,NULL}, idx_(idx)
{
    if( parent_ != NULL )
        this->level_ = parent_->GetLevel() + 1;
    else
        this->level_ = 0;
}

/** @brief Manages the destruction of Node objects and memory management.
 * @details Moooar
 */
Node::~Node() {
    for( size_t ic(0); ic<4; ic++) {
        if( this->children_[ic] != NULL ) {
            delete this->children_[ic];
        }
    }
}

/** @brief Returns the parent node of the current object.
 */
Node* Node::GetParent() const {
    return parent_;
}

/** @brief Returns the child node stored at the input index.
 * @param[in] child_idx Index of the child node in the current node's array.
 */
Node* Node::GetChild( short int child_idx ) const {
    return children_[child_idx];
}

/** @brief Finds the non-empty Node on the same level.
 * @details During the domain decomposition, initiated by the
 * RootNode::Decomposition method, we recursively generated quad-tree decomposed
 * subdomains, each associated with a Node located at a lower level of
 * refinement. @n
 * Some procedures within the SURMISE code require browsing all the Node objects
 * from a specific level (for example, when building the first multipoles at the
 * LeafNode level). A method navigating through all these nodes sequentially is
 * thus needed, this is the aim of the current function. @n
 * When a Node calls GetNext, the parent Node is inspected, and children with a
 * higher index are looked for and returned. If the current node has index 3,
 * no antecedant node on the same level, and contained in the same quad-tree
 * subdomain exist. A recursive search is then launched and navigates up the
 * tree diagram until a larger, neighboring subdomain with greater index than
 * the current one is found. This subdomain is entered and explored until a Node
 * matching the caller's level is found, which the method returns. @n
 * Node::GetNext returns NULL when the recursion tries to look for the
 * RootNode's parent, indicating no Node with greater index can be found.
 * @returns A pointer to the next Node found on the same level.
 */
Node* Node::GetNext() const {
    // If we're at the Root, we want to return immediately.
    if( this->level_ == 0 )
        return NULL;
    // Node is not last in quad tree decomposition. Parent has a child with
    // higher index (which we return, most common case).
    unsigned int idx(this->index_%4);
    if( idx<3 )
        return this->GetParent()->GetChild( idx + 1 );
    // Else, we need to search the next Node among the parents, get the current
    // node's.
    Node* nxt( this->GetParent() );
    // We stepped a level higher. Make sure we return the next node on the same
    // level.
    while( true ) {
        // If reached RootNode, return NULL; we come from the last child of this
        // RootNode and cannot look further.
        if( nxt->GetLevel() == 0)
            return NULL;
        // This node's parent contains more interesting children, select the
        // next one (horizontal move in the tree).
        if( nxt->GetIndex() < 3) {
            nxt = nxt->GetParent()->GetChild( nxt->GetIndex() + 1 );
            // We are in a new branch, follow the first child until we're back
            // at the starting level.
            while( nxt->GetLevel() != this->level_ )
                nxt = nxt->GetChild( 0 );
            return nxt;
        }
        // If this statement is reached, this means the investigated parent is
        // also the last child of the "upper" level, we need to search higher up
        // the tree.
        nxt = nxt->GetParent();
    }
}

/** @brief Returns the edge sharing Node to the caller in the requested
 * direction.
 * @details When discretizing the domain, a recursive quad-tree is built. This
 * approach leads to separation between nodes in memory to be very different
 * from their physical separation in the 2D plane in some cases. @n
 * Consider a Node which's corner is at the center of the simulated region. The
 * three other nodes sharing this corner lay in well separated regions of the
 * quad-tree, since the first domain separation that was performed split these
 * four regions. Hence, the caller node's neighbour would only be found by
 * browsing the tree upwards until the RootNode, and searching downwards after
 * that. @n
 * Clearly, bruteforcing this search is not a good idea, hence the current
 * function. It takes a direction and returns a pointer to the node neighbouring
 * the caller in this direction. In the event that the node cannot be found
 * (encounters the simulated domain's border), a NULL pointer is returned.
 * @param[in] direction Enumerated type describing all possible directions in a
 * two dimensional plane, up/down/left/right.
 * @returns A pointer to the nearest node in the specified direction from the
 * calling node.
 */
Node* Node::Move( ZDIR direction ) const {
    // If at the Root, can't move, return NULL
    if( this->level_ == 0 )
        return NULL;
    // Compute the targeted index
    int closest_idx( (int)this->idx_ + direction );
    int target_idx(0);
    if( direction % 2 == 0 ) {
        // The direction is LEFT/RIGHT
        // Required conversion (0,1)<->(2,3), ie flipping the 2nd bit
        target_idx = this->idx_ ^ 1UL << 1;
    } else {
        // The direction is UP/DOWN
        // Required conversion (0,2)<->(1,3), ie flipping the 1st bit
        target_idx = this->idx_ ^ 1UL << 0;
    }
    if( closest_idx == target_idx ) {
        // The target index is within the current parent's reach. Can return
        // directly
        return this->GetParent()->GetChild( target_idx );
    }
    // Else, the requested node is separated by a border defined higher up the
    // tree. Need to perform the index conversions (0,2)<->(1,3) for up/down
    // moves and (0,1)<->(2,3) for left/right ones. At first, one could think a
    // modulo 4 would have been ssufficient, but problems are appearing (for
    // example, start idx = 3, move up, so mod 4 says zero...)
    // Try the same move on the parent node.
    Node* neighbour( this->GetParent()->Move(direction) );
    // If we reached the Root to realize we're against a wall, return NULL
    if( neighbour == NULL )
        return NULL;
    return neighbour->GetChild( target_idx );
}

/** Is nearest neighbour?
 * @todo dx definition dirty.
 */
bool Node::IsNN( const Node* other ) const {
    std::array<double,4> obnds( other->GetBounds() );
    double dx( this->xybnds[1] - this->xybnds[0] + 1.e-6);
    if( std::abs( this->xybnds[0] - obnds[0] ) < dx
            and std::abs( this->xybnds[2] - obnds[2] ) < dx
            and this!=other)
        return true;
    return false;
}

/** @brief Build an array of pointers to the caller's nearest neighbours.
 * @details This method returns a fixed-size array of pointers to nodes
 * adjacent (edge and corner sharing) to the calling node. These nodes are
 * determined using the Node::Move method, so that when encountering walls (i.e.
 * no nodes), a NULL pointer is returned.
 * @returns A fixed-size array containing the caller node's nearest neighbours.
 */
std::array<Node*,8> Node::GetNN() const {
    // Take ordering convention
    // 2  4   7
    // 1 this 6
    // 0  3   5
    std::array<Node*,8> nn;
    nn.fill( NULL );
    nn[1] = this->Move(LEFT);
    if(nn[1] != NULL) {
        nn[0] = nn[1]->Move(DOWN);
        nn[2] = nn[1]->Move(UP);
    }
    nn[3] = this->Move(DOWN);
    nn[4] = this->Move(UP);
    nn[6] = this->Move(RIGHT);
    if( nn[6] != NULL ) {
        nn[5] = nn[6]->Move(DOWN);
        nn[7] = nn[6]->Move(UP);
    }
    return nn;
}

/** @brief Searches and returns the current node's interaction list.
 * @details In the FMM, the multipole to local translation implies the
 * propagation of a node's mutlipole to its interaction list, that is, to all
 * the child nodes located within the children of its parent nearest neighbours,
 * while excluding the nodes that are directly then nearest neighbours of the
 * caller. This list amounts to a maximum of 27 elements. The interaction list
 * array is first filled with NULL pointers, which are replaced as members of
 * the interaction list are found.
 * @returns A fixed size array of nodes in the caller's interaction list.
 */
std::array<Node*,27> Node::GetIN() const {
    std::array<Node*,8> parentNN( this->GetParent()->GetNN() );
    std::array<Node*,27> in;
    in.fill( NULL );
    unsigned int ii(0);
    Node* child(NULL);
    for( unsigned int ip(0); ip<8; ip++ ) {
        if( parentNN[ip] != NULL ) {
            for( unsigned int ic(0); ic<4; ic++ ) {
                child = parentNN[ip]->GetChild( ic );
                if( not this->IsNN(child) ) {
                    in[ii] = child;
                    ii++;
                }
            }
        }
    }
    return in;
}

/** @brief Instanciates a child to the current Node.
 * @details This methods dynamically allocates a new node at position `child`
 * in the Node's children_ array.
 * @warning This method is not memory safe. Any already allocated child node is
 * not freed.
 * @todo Obviously check memory if performance not at stake.
 * @param[in] child Index of the child Node to instanciate (0-3).
 * @returns An error code.
 */
SError Node::InitChild( int child ) {
    this->children_[child] = new Node(this,child);
    bool low_quad(child==0 or child==2), left_quad(child==0 or child==1);
    this->children_[child]->SetBounds(
            this->xybnds[0] + 0.5 * (double)(not left_quad) * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[1] - 0.5 * (double)left_quad * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[2] + 0.5 * (double)(not low_quad) * ( this->xybnds[3] - this->xybnds[2] ),
            this->xybnds[3] - 0.5 * (double)low_quad * ( this->xybnds[3] - this->xybnds[2] )
            );
    return E_SUCCESS;
}

/** @brief Instanciates all children to the current Node.
 * @details This methods dynamically allocates four new nodes in the Node's
 * children_ array.
 * @warning This method is not memory safe. Any already allocated child node is
 * not freed.
 * @todo Obviously check memory if performance not a stake.
 * @returns An error code.
 */
SError Node::InitAllChildren() {
    for( int ic(0); ic<4; ic++ ) {
        this->InitChild( ic );
    }
    return E_SUCCESS;
}

/** @brief Instanciates a child to the current Node.
 * @details This methods dynamically allocates a new node at position `child`
 * in the Node's children_ array.
 * @warning This method is not memory safe. Any already allocated child node is
 * not freed.
 * @todo Obviously check memory if performance not at stake.
 * @param[in] child Index of the child Node to instanciate (0-3).
 * @returns An error code.
 */
SError Node::InitLeaf( int child ) {
    this->children_[child] = new LeafNode(this,child);
    bool low_quad(child==0 or child==2), left_quad(child==0 or child==1);
    this->children_[child]->SetBounds(
            this->xybnds[0] + 0.5 * (double)(not left_quad) * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[1] - 0.5 * (double)left_quad * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[2] + 0.5 * (double)(not low_quad) * ( this->xybnds[3] - this->xybnds[2] ),
            this->xybnds[3] - 0.5 * (double)low_quad * ( this->xybnds[3] - this->xybnds[2] )
            );
    return E_SUCCESS;
}

/** @brief Instanciates all children to the current Node.
 * @details This methods dynamically allocates four new nodes in the Node's
 * children_ array.
 * @warning This method is not memory safe. Any already allocated child node is
 * not freed.
 * @todo Obviously check memory if performance not a stake.
 * @returns An error code.
 */
SError Node::InitAllLeafs() {
    for( int ic(0); ic<4; ic++ ) {
        this->InitLeaf( ic );
    }
    return E_SUCCESS;
}

/** @brief Computes the geometrical center of the calling node from its
 * boundaries.
 */
std::array<double,2> Node::Center() const {
    return this->com;
    /*return std::array<double,2>{(double)0.5*(this->xybnds[0]+this->xybnds[1]),
    (double)0.5*(this->xybnds[2]+this->xybnds[3])};*/
}

/** @brief Returns an array containing the node's geometrical bounds.
 */
std::array<double,4> Node::GetBounds() const {
    return std::array<double,4>{
        this->xybnds[0],
        this->xybnds[1],
        this->xybnds[2],
        this->xybnds[3]};
}

/** @brief Sets the (physical) boundaries of the domain
 * @param[in] x1,x2,y1,y2 Boundaries of the domain
 * @returns An error code
 */
SError Node::SetBounds( double x1, double x2, double y1, double y2 ) {
    this->xybnds[0] = x1;
    this->xybnds[1] = x2;
    this->xybnds[2] = y1;
    this->xybnds[3] = y2;
    return E_SUCCESS;
}

/** @brief Splits the current domain into four parts.
 * @details The splitting is performed recursively until enough levels have been
 * created. Pointers to Particle objects are passed from the RootNode
 * recursively to until a LeafNode is generated. These pointers are not sotred
 * within the class by this method.
 * @param[in] parts Pointers to particles located in the current domain. They
 * will be separated in arrays corresponding to the generated subdomains and
 * sent down the tree by recursively calling Node::Decompose.
 * @returns Integer containing an error code. If the method exits cleanfully,
 * zero is returned.
 */
SError Node::Decompose( std::vector<Particle*> &parts ) {
    // Create children and further decompose
    // If required accuracy is obtained, create LeafNode children
    if(this->level_ < this->maximum_level) {
        this->InitAllChildren();
    } else {
        this->InitAllLeafs();
    }
    for( int ic(0); ic<4; ic++ ) {
        std::vector<Particle*> subparts;
        for( unsigned int ip(0); ip<parts.size(); ip++ ){
            if( this->GetChild(ic)->BelongsTo(parts[ip]) ) {
                subparts.push_back(parts[ip]);
            }
        }
        this->GetChild(ic)->Decompose(subparts);
    }
    return E_SUCCESS;
}

/** @brief Performs a time step for the temporal evolution.
 * @details This function is called recursively from the RootNode and expects
 * the LeafNode objects to contains relevant masses (i.e., the domain has been
 * decomposed already). It also assumes relevant forces have been computed by
 * the Fast Multipole Algorithm.
 * @param[in] dt Time step of the temporal evolution.
 * @warning Choosing a high time step will result in unreliable results, or to
 * unstable simulations. The temporal evolution algorithm being *explicit*, a
 * high `dt` is likely to be problematic.
 * @returns Integer containing an error code. If the method exits cleanfully,
 * zero is returned.
 */
SError Node::TimeEvolution( double dt ) {
    if( this->GetLevel() > 2) {
        std::array<double,2> upstream_force(this->GetParent()->force);
        this->force[0] += upstream_force[0];
        this->force[1] += upstream_force[1];
    }
    for( unsigned int ic(0); ic<4; ic++ ) {
        if( this->children_[ic] != NULL )
            this->children_[ic]->TimeEvolution( dt );
    }
    return E_SUCCESS;
}

/** Recursively calls leaf nodes' implementation of the method to update the
 * positions and velocities of the simulated particles.
 * @param[in] dt Time step for the temporal evolution.
 * @returns An error code
 */
SError Node::TimeEvolutionMasses(double dt) {
    for( auto c : this->children_ )
        c->TimeEvolutionMasses(dt);
    return E_SUCCESS;
}

/** @brief  Reassigns a Particle to a new subdomain.
 * @details After a time evolution step is performed, the particles may hop out
 * of the LeafNode's domains they were assigned to. In order to keep the FMM
 * algorithm running correctly, those Particle(s) need to be reassigned to the
 * correct leaf. The leaf nodes will routinely call this method using a virtual
 * redirection. @n
 * When receiving a node pointer, the method checks if the particle is within
 * its covered area. If not, it passes it to its parent node. Otherwise, it
 * iterates over its children nodes, identifies the one hosting the particle,
 * and does so until a leaf is reached, and the particle assigned.
 * @todo Improve the virtual redirection of the leaf, or even remove it: as long
 * as the first llines of Node::Reassign still checks if invokking the parent is
 * needed, we can avoid an override.
 * @returns An error code.
 */
SError Node::Reassign( Particle* p ) {
    // If the particle is not in this domain, must pass to parent and enlarge
    // search.
    if( not this->BelongsTo(p) ) {
        this->GetParent()->Reassign(p);
    } else {
        // The particle is in this domain, find which child should get it.
        Node* recv_node(this);
        while( recv_node->GetLevel() != maximum_level+1 ) {
            for( unsigned int ic(0); ic<4; ic++ ) {
                if( recv_node->GetChild(ic) != NULL ) {
                    if( recv_node->GetChild(ic)->BelongsTo(p) ) {
                        // Found the child, shift pointer, break the current
                        // level's loop and try again until leaf level.
                        recv_node = recv_node->GetChild(ic);
                        break;
                    }
                }
            }
        }
        // We're at the last level, add the particle to the leaf's stack.
        recv_node->AddParticle(p);
    }
    return E_SUCCESS;
}

/** @brief Indirection to the leaf's method.
 * @note Since a node cannot store particle pointers, we assume a reassignement
 * was meant and call the method from the current Node.
 * @returns An error code.
 */
SError Node::AddParticle( Particle* p ) {
    // A node cannot store particles, redirects to reassignment method.
    this->Reassign(p);
    return E_SUCCESS;
}

/** @brief Computes the force between two given nodes.
 * @details This method updates the current node's forces by adding the
 * contribution of a distant node (passed as parameter). The force is computed
 * between the nodes' _geometrical_ centers. It is based only on the total mass
 * of the nodes (and not a multipole expansion of these). In principle, nodes
 * from different levels may call this method.
 * @note Although gravitational interaction is a symmetric process, this method
 * only updates the force of the caller node. A separate call is needed if the
 * remote node requires the influence of the caller node to be accounted for.
 * @param[in] other Pointer to a remote node which's interaction must be
 * accounted for.
 * @todo Try to symmetrize force updates and consider their center of mass
 * instead of their geometrical one.
 * @returns An error code.
 */
SError Node::Interact( const Node* other ) {
    std::array<double,2> tcen(this->Center());
    std::array<double,2> ocen(other->Center());
    double dist2( (tcen[0]-ocen[0])*(tcen[0]-ocen[0]) + (tcen[1]-ocen[1])*(tcen[1]-ocen[1]) );
    this->force[0] += - other->mass / dist2 * (tcen[0]-ocen[0]);
    this->force[1] += - other->mass / dist2 * (tcen[1]-ocen[1]);
    return E_SUCCESS;
}

/** @brief Updates the total mass of the current Node.
 * @details Computes the total mass of the caller Node by summing the mass of
 * its children. This function calls the method Node::GetMass from its children,
 * thus, the mass of the children is not updated nor modified by this call. The
 * children nodes thus need to be updated for the opreations performed by this
 * function to make sense. This is done by calling Node::GatherMasses for all
 * children prior to calling this function.
 * @returns An error code.
 */
SError Node::GatherMasses() {
    this->mass = this->children_[0]->GetMass() + this->children_[1]->GetMass() +this->children_[2]->GetMass() + this->children_[3]->GetMass();
    if( this->mass != 0.0 ) {
        this->com[0] = ( this->children_[0]->GetMass() * this->children_[0]->com[0] + this->children_[1]->GetMass() * this->children_[1]->com[0] + this->children_[2]->GetMass() * this->children_[2]->com[0] + this->children_[3]->GetMass() * this->children_[3]->com[0] ) / this->mass;
        this->com[1] = ( this->children_[0]->GetMass() * this->children_[0]->com[1] + this->children_[1]->GetMass() * this->children_[1]->com[1] + this->children_[2]->GetMass() * this->children_[2]->com[1] + this->children_[3]->GetMass() * this->children_[3]->com[1] ) / this->mass;
    } else {
        this->com = std::array<double,2>{(double)0.5*(this->xybnds[0]+this->xybnds[1]),
            (double)0.5*(this->xybnds[2]+this->xybnds[3])};
    }
    return E_SUCCESS;
}

/** Accessor to the current node's mass.
 * @returns The node's total mass.
 */
double Node::GetMass() const { return mass; }

/** Sets a node's total mass. This methods just sets the mass and does not
 * perform any other kind of computation.
 * @param[in] m Mass that needs to be registered in the node object.
 * @returns An error code.
 */
SError Node::SetMass( double m ){
    this->mass = m;
    return E_SUCCESS;
}

/** Sets the caller node's forces to zero. This method is called for the whole
 * tree before starting time evolution computations.
 * @returns An error code.
 */
SError Node::ResetForces() {
    this->force[0] = 0.;
    this->force[1] = 0.;
    return E_SUCCESS;
}

/** Adds a 2D force to the current node.
 * @param[in] upstream_force Reference to a size-2 array.
 * @returns An error code.
 */
SError Node::AddForce( const std::array<double,2>& upstream_force ) {
    this->force[0] += upstream_force[0];
    this->force[1] += upstream_force[1];
    return E_SUCCESS;
}

/** Accessor to the caller node's force.
 * @todo If performance hit, remove the array copy for something fancier.
 * @returns A copy of the force acting on the node's center.
 */
std::array<double,2> Node::GetForce() const {
    return this->force;
}

SError Node::SetCOM( const std::array<double,2>& cen ) {
    this->com = cen;
    return E_SUCCESS;
}

/** Returns the level of the current Node. The RootNode (likely a parent of the
 * current one) has level zero.
 */
unsigned int Node::GetLevel() const { return level_; }

/** Returns the current Node's index in its local sub-tree.
 * @returns The Node's index.
 */
unsigned int Node::GetIndex() const { return idx_; }

/** @brief Returns true if the argmuent Particle is in the domain
 * @details Many more things
 * @param[in] p Pointer to a Particle object.
 * @returns Boolean indicating the Particle's belonging to the domain.
 */
bool Node::BelongsTo( Particle *p ) const {
    double xp(p->pos[0]); double yp(p->pos[1]);
    if( xp > this->xybnds[0] && xp < this->xybnds[1] &&
            yp > this->xybnds[2] && yp < this->xybnds[3] )
        return true;
    return false;
}

//////////////////// ROOTNODE ///////////////////////
/** @brief Creates the simulation's root node
 * @details This call create the first node of the simulation, that is,
 * instantiates the whole smiulation domain/box. The Particle(s) are stored in
 * the configuration parameter.
 * @param[in] cf Input configuration for the simulation. Most parameters it
 * contains will be "unpacked" in the class' members.
 */
RootNode::RootNode( SConfig& cf )
    :Node(NULL,0), conf_(cf)
{
    this->SetBounds(0.,this->conf_.dsize,0.,this->conf_.dsize);
}

/** @brief Launches the requested simulation.
 * @details After the initialization of the RootNode class, a call to this
 * method will launch a simulation with criteria matching those imposed in the
 * configuration object RootNode::conf_. The simulation is just a series of
 * step, with a step looking like:
 *
 * 1. Domain decomposition and Particle objects assignements,
 * 2. Fast Multipole Method for force computation,
 * 3. Time evolution for each Particle,
 * 4. [Optional] Upwards update to the RootNode,
 * 5. I/O Routines,
 * 6. Verify end of run conditions.
 *
 * Once the number of demanded iterations has been reached, or once the
 * allocated time for the simulation is reached, a final I/O operation is
 * performed to generate restart files, and the function exits with code zero.
 *
 * @returns Integer containing an error code. If the method exits cleanfully,
 * zero is returned.
 * @todo Verify if step 4 in the step description is needed.
 */
SError RootNode::Run() {
    SError run_err(E_SUCCESS);
    this->Decompose( this->conf_.parts );
    int max_it = this->conf_.max_iter;
    for( int it(0); it<max_it; it++ ) {
        // Particles are distributed, perform FMM
        run_err = this->TimeEvolution( this->conf_.dt );
        run_err = this->TimeEvolutionMasses( this->conf_.dt );
        /// @todo Error management, check computational cost.
        for(auto p : this->conf_.parts ) {
            std::cout << std::setprecision(10) << p->id << "," << p->mass << "," << p->pos[0] << "," << p->pos[1] << "," << p->vel[0] << "," << p->vel[1] << "," << p->frc[0] << "," << p->frc[1] << std::endl;
        }
    }
    return run_err;
}

/** @brief Initiates the simulation's domain tree decomposition.
 * @details Computes the required levels of the tree based on the number of
 * particles, and recursively distributes Particle objects to the child nodes.
 * @todo Update mechanism to handle adaptative structure
 * @param[in] parts Particle(s) that  are to be distributed in the tree, most
 * likely
 * those contained in the simulation configuration, RootNode::conf_.
 * @returns Integer containing an error code. If the method exits cleanfully,
 * zero is returned.
 */
SError RootNode::Decompose( std::vector<Particle*> &parts ) {
    // Compute the desired amount of refinement
    Node::maximum_level = (int) std::max(
            2.,
            ///@todo Check whether -1 needed here
            std::ceil(std::log( this->conf_.npart )) - 1
            );
    // Create children and further decompose
    this->InitAllChildren();
    for( int ic(0); ic<4; ic++ ) {
        std::vector<Particle*> subparts;
        for( unsigned int ip(0); ip<this->conf_.npart; ip++ ){
            if( this->GetChild(ic)->BelongsTo(this->conf_.parts[ip]) ) {
                subparts.push_back(this->conf_.parts[ip]);
            }
        }
        this->GetChild(ic)->Decompose(subparts);
    }
    return E_SUCCESS;
}

/** @brief Initiates a temporal evolution of a given time step.
 * @details This method calls the class' children nodes time evolution
 * mechanism. It assumes Fast Multipole Aprroximation of forces have been
 * performed previously.
 * @param[in] dt Time step of the temporal evolution.
 * @warning Choosing a high time step will result in unreliable results, or to
 * unstable simulations. The temporal evolution algorithm being *explicit*, a
 * high `dt` is likely to be problematic.
 * @returns Integer containing an error code. If the method exits cleanfully,
 * zero is returned.
 */
SError RootNode::TimeEvolution( double dt ) {
    // Upwards sweep, gather masses
    for(int iL(this->maximum_level+1); iL>1; iL--) {
        Node* ptr(this);
        // Move pointer down the tree (to the expected level)
        while(ptr->GetLevel()!=iL){
            ptr = ptr->GetChild(0);
        }
        // Iterate on all level nodes
        do {
            // Gather children masses
            ptr->ResetForces();
            ptr->GatherMasses();
            ptr = ptr->GetNext();
        }while( ptr != NULL );
    }
    // Downwards sweep, compute interaction on each level
    for(int iL(2); iL<this->maximum_level+2; iL++) {
        Node* ptr(this);
        // Move pointer down the tree (to the expected level)
        while(ptr->GetLevel()!=iL){
            ptr = ptr->GetChild(0);
        }
        do{
            // Compute interaction of each cell with far NN.
            std::array<Node*,27> fnn(ptr->GetIN());
            for( auto dnode : fnn ) {
                // Add the interaction of the neighbour if it exists
                if(dnode!=NULL) {
                    ptr->Interact(dnode);
                } else {
                    // If we reached a NULL, the rest of the table is empty.
                    break;
                }
            }
            ptr = ptr->GetNext();
        }while(ptr!=NULL);
    }
    // Call children time evolution.
    for( unsigned int ic(0); ic<4; ic++)
        if( this->GetChild(ic) != NULL )
            this->GetChild(ic)->TimeEvolution(dt);
    return E_SUCCESS;
}

/** @brief Reassignement procedure for the most shallow level.
 * @details If the RootNode is asked to reassign a Particle, it checks whether
 * it is located in the simulated area at all. If it finds a child hosting the
 * particle, it is passed down again, see Node::Reassign. Else, the method ends
 * _without reassigning the particle_.
 * @note Particle objects hopping out of the simulated area are not reassigned,
 * and considered lost. May be revised in future versions.
 * @returns An error code.
 */
SError RootNode::Reassign( Particle* p ) {
    for( unsigned int ic(0); ic<4; ic++ ) {
        if( this->GetChild(ic)->BelongsTo(p) ) {
            // Found the expected child, send the particle there.
            this->GetChild(ic)->Reassign(p);
            break;
        }
    }
    // XXX - Reaching this statement, the particle is possibly not assigned to
    // any domain (i.e. outside the simulated area), and remains unhandled for
    // the rest of the simulation (it was deleted from its leaf when
    // reassignment began).
    return E_SUCCESS;
}

//////////////////// LEAFNODE ///////////////////////
/** @brief Builds a node loacted at the bottom of the tree, i.e. a leaf.
 * @details The node is build in the usual way, except we require a Particle
 * array to be passed. It is not possible to modify the LeafNode::children_
 * member, that will be left empty.
 * @param[in] parent Pointer to a parent node.
 * @param[in] child Index of the child within its parent's children subtree.
 */
LeafNode::LeafNode( Node *parent, unsigned int child )
    :Node(parent,child)
{
    // Make sure the storage is not initialized.
    this->subparts_.empty();
}

/** @brief Builds a node loacted at the bottom of the tree, i.e. a leaf.
 * @details The node is build in the usual way, except we require a Particle
 * array to be passed. It is not possible to modify the LeafNode::children_
 * member, that will be left empty.
 * @param[in] parent Pointer to a parent node.
 * @param[in] child Index of the child within its parent's children subtree.
 * @param[in] parts Array of Particle objects that are located within the node's
 * boundaries.
 */
LeafNode::LeafNode( Node *parent, unsigned int child, std::vector<Particle*> &parts )
    :Node(parent,child), subparts_(parts)
{}

/** @brief Ends the decomposition procedure.
 * @details When recursive calls to Node::Decompose have been performed, a
 * satisfactory level of decomposition is eventually achieved. When this
 * happens, Node::Decompose will generate the last level with leaf nodes and
 * call their decomposition method. It simply stores the assigned particle in
 * the subparts_ array.
 * @param[in] parts Particle(s) that are to be distributed in the tree, most likely
 * those contained in the simulation configuration, RootNode::conf_.
 * @returns Integer containing an error code. If the method exits cleanfully,
 * zero is returned.
 */
SError LeafNode::Decompose( std::vector<Particle*> &parts ) {
    this->subparts_ = parts;
    return E_SUCCESS;
}

/** @brief Handles the system's time evolution.
 * @details This method is called after recursive calls have been made to
 * Node::TimeEvolution from higher levels. An explicit scheme is applied to each
 * mass and evolves their velocity and position given the force acting on them.
 * @note Before calling this method, the forces acting upon the particles should
 * be computed using the FMM.
 * @todo Why iteratively moving to lower levels for time evolution? As long as
 * the memory is local, this task may only be performed by the RootNode.
 * @param[in] dt Time step of the temporal evolution.
 * @warning Choosing a high time step will result in unreliable results, or to
 * unstable simulations. The temporal evolution algorithm being *explicit*, a
 * high `dt` is likely to be problematic.
 * @returns Integer containing an error code. If the method exits cleanfully,
 * zero is returned.
 */
SError LeafNode::TimeEvolution( double dt ) {
    if( this->GetLevel() > 2) {
        std::array<double,2> upstream_force(this->GetParent()->GetForce());
        this->AddForce( upstream_force );
    }
    for( unsigned int ip(0); ip<this->subparts_.size(); ip++ ) {
        Particle* p(this->subparts_[ip]);
        // Add upstream force (Barnes-Hut)
        std::array<double,2> mfrc(this->GetForce());
        p->frc[0] = p->mass*mfrc[0]; p->frc[1] = p->mass*mfrc[1];
        // Add force from particles within same region
        for( unsigned int ip2(0); ip2<this->subparts_.size(); ip2++ ) {
            Particle* other(this->subparts_[ip2]);
            if( p!=other ) {
                std::array<double,2> ppfrc(p->PForce(other));
                p->frc[0] += ppfrc[0];
                p->frc[1] += ppfrc[1];
            }
        }
        // Add force from nearest neighbours
        std::array<Node*,8> nn(this->GetNN());
        for( int inn(0); inn<8; inn++ ){
            if(nn[inn]!=NULL) {
                std::vector<Particle*> nnp(nn[inn]->GetParticles());
                for( unsigned int ip2(0); ip2<nnp.size(); ip2++ ){
                    if(p!=nnp[ip2]){
                        std::array<double,2> ppfrc(p->PForce(nnp[ip2]));
                        p->frc[0] += ppfrc[0];
                        p->frc[1] += ppfrc[1];
                    }
                }
            }
        }
    }
    return E_SUCCESS;
}

/** @brief Updates a leaf's particles positions and velocities.
 * @details This method applies the forces computed by LeafNode::TimeEvolution
 * to the LeafNode's particles for a timestep 'dt'. Once the update is
 * performed, the particles' positions are checked to see whether a domain
 * reassignement is required. Such a function allows decoupling the force
 * computation and the particle's temporal evolution.
 * @param[in] dt Time step for the temporal evolution.
 * @returns An error code.
 */
SError LeafNode::TimeEvolutionMasses(double dt) {
    for( unsigned int ip(0); ip<this->subparts_.size(); ip++ ) {
        Particle* p(this->subparts_[ip]);
        p->vel[0] += p->frc[0] / p->mass * dt;
        p->vel[1] += p->frc[1] / p->mass * dt;
        p->pos[0] += p->vel[0] * dt;
        p->pos[1] += p->vel[1] * dt;
        if( not this->BelongsTo(p) ) {
            // Remove the reference to the particle from the leaf and reassign.
            this->subparts_.erase( this->subparts_.begin() + ip );
            this->Reassign(p);
        }
    }
    return E_SUCCESS;
}

/** @brief Reassigns a Particle to a new domain. The parent function is called.
 * @param[in] p A pointer to a (local) Particle.
 * @returns An error code indicating if the method succeeded.
 */
SError LeafNode::Reassign( Particle* p ) {
    this->GetParent()->Reassign( p );
    return E_SUCCESS;
}

/** @brief Adds the particle to the leaf's storage
 * @note Performs a simple pushback. We assume an antecedent check was performed
 * to verify the Particle is in the simulated area.
 * @returns An error code.
 */
SError LeafNode::AddParticle( Particle* p ) {
    // Pushback without asking questions.
    this->subparts_.push_back(p);
    return E_SUCCESS;
}

/** Sets the leaf's mass as the sum of the particles it contains.
 * @returns An error code.
 */
SError LeafNode::GatherMasses() {
    this->SetMass(0.);
    double m(0.);
    for( auto p : subparts_ )
        m+= p->mass;
    this->SetMass( m );
    std::array<double,2> cen({0.,0.});
    if(this->GetMass()==0.0){
        std::array<double,4> xybnds = this->GetBounds();
        cen = std::array<double,2>{(double)0.5*(xybnds[0]+xybnds[1]),
            (double)0.5*(xybnds[2]+xybnds[3])};
    } else {
        for( auto p : subparts_ ) {
            cen[0] += p->mass * p->pos[0] / this->GetMass();
            cen[1] += p->mass * p->pos[1] / this->GetMass();
        }
    }
    this->SetCOM(cen);
    return E_SUCCESS;
}

/** @brief Returns the current leaf's particles list.
 * @details This method aims to provide access to particles living in a given
 * leaf. It is useful when computing the nearest neighbor forces, as in
 * LeafNode::TimeEvolution.
 * @todo Remainder to check whether this is a performance issue. If yes, provide
 * "sequential" accessor to avoid returning the whole particle container (or
 * return a reference, or whatever, this copy is just plain dirty).
 * @returns A copy of the pointers to particles stored in the caller leaf.
 */
std::vector<Particle*> LeafNode::GetParticles() const {
    return this->subparts_;
}
