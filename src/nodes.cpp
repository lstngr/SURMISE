#include <iostream>
#include <cmath>
#include "nodes.hpp"

// Declare the static member here (not in the header!)
unsigned int Node::maximum_level = 2;

//////////////////// NODE ///////////////////////

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
    if( this->idx_<3 )
        return this->GetParent()->GetChild( this->idx_ + 1 );
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
    std::array<float,4> obnds( other->GetBounds() );
    float dx( this->xybnds[1] - this->xybnds[0] + 1.e-6);
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
            this->xybnds[0] + 0.5 * (float)(not left_quad) * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[1] - 0.5 * (float)left_quad * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[2] + 0.5 * (float)(not low_quad) * ( this->xybnds[3] - this->xybnds[2] ),
            this->xybnds[3] - 0.5 * (float)low_quad * ( this->xybnds[3] - this->xybnds[2] )
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
            this->xybnds[0] + 0.5 * (float)(not left_quad) * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[1] - 0.5 * (float)left_quad * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[2] + 0.5 * (float)(not low_quad) * ( this->xybnds[3] - this->xybnds[2] ),
            this->xybnds[3] - 0.5 * (float)low_quad * ( this->xybnds[3] - this->xybnds[2] )
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
std::array<float,2> Node::Center() const {
    return std::array<float,2>{(float)0.5*(this->xybnds[0]+this->xybnds[1]),
    (float)0.5*(this->xybnds[2]+this->xybnds[3])};
}

/** @brief Returns an array containing the node's geometrical bounds.
 */
std::array<float,4> Node::GetBounds() const {
    return std::array<float,4>{
        this->xybnds[0],
        this->xybnds[1],
        this->xybnds[2],
        this->xybnds[3]};
}

/** @brief Sets the (physical) boundaries of the domain
 * @param[in] x1,x2,y1,y2 Boundaries of the domain
 * @returns An error code
 */
SError Node::SetBounds( float x1, float x2, float y1, float y2 ) {
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
    ///@todo Implement time evolution algorithm.
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
    float xp(p->pos[0]); float yp(p->pos[1]);
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
    // TEST ZONE //
    Node* lowest(this);
    while( lowest->GetLevel() != this->maximum_level+1 )
        lowest = lowest->GetChild(0);
    while( lowest!=NULL ) {
        std::array<float,4> bnds(lowest->GetBounds());
        Node* upper( lowest->Move(RIGHT) );
        int idxupper(-1);
        if( upper != NULL ) {
            idxupper = upper->GetIndex();
        }
        // std::cout << "[LEAF] Level(" << lowest->GetLevel() << ") - Index(" << lowest->GetIndex() << ") - Center(" << 0.5*(bnds[0]+bnds[1]) << "," << 0.5*(bnds[2]+bnds[3]) << ") - Right neighbour - (" << idxupper << ")" << std::endl;
        std::array<Node*,27> smthg(lowest->GetIN());
        int ii(0);
        while( smthg[ii]!=NULL && ii<27 ) {
            std::array<float,2> cen(smthg[ii]->Center());
            //std::cout << "Center (" << cen[0] << "," << cen[1] << ") ";
            ii++;
        }
        //std::cout << std::endl;
        lowest = lowest->GetNext();
    }
    ///////////////
    for( int it(0); it<max_it; it++ ) {
        // Particles are distributed, perform FMM
        run_err = this->TimeEvolution( this->conf_.dt );
        /// @todo Error management, check computational cost.
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
    return E_SUCCESS;
}
