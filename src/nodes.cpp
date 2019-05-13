#include <cmath>
#include "nodes.hpp"

// Declare the static member here (not in the header!)
int Node::maximum_level = 2;

//////////////////// NODE ///////////////////////

/** @brief Creates a generic node.
 * @param[in] parent Pointer to a parent node. If none is provided, the pointer
 * is left null and the node is interpreted to be at the root.
 */
Node::Node( Node *parent )
    :parent_(parent), children_{NULL,NULL,NULL,NULL}
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

/** @brief Instanciates a child to the current Node.
 * @details This methods dynamically allocates a new node at position `child`
 * in the Node's children_ array.
 * @warning This method is not memory safe. Any already allocated child node is
 * not freed.
 * @todo Obviously check memory if performance not a stake.
 * @param[in] child Index of the child Node to instanciate (0-3).
 * @returns An error code.
 */
SError Node::InitChild( int child) {
    this->children_[child] = new Node(this);
    ///@todo Implement correct bounds
    this->children_[child]->SetBounds(
            this->xybnds[0] + 0.5 * (child==0 or child==2) * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[1] - 0.5 * (child==1 or child==3) * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[2] + 0.5 * (child==0 or child==1) * ( this->xybnds[1] - this->xybnds[0] ),
            this->xybnds[3] - 0.5 * (child==2 or child==3) * ( this->xybnds[1] - this->xybnds[0] )
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
    ///@todo Implement tree decomposition.
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
    :Node(NULL), conf_(cf)
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
            std::ceil(std::log( this->conf_.npart ))
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
        this->GetChild(ic)->Decompose(parts);
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
 * @param[in] parts Array of Particle objects that are located within the node's
 * boundaries.
 */
LeafNode::LeafNode( Node *parent, std::vector<Particle*> parts )
    :Node(parent), subparts_(parts)
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
