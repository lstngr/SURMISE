#include <iostream>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "nodes.hpp"

//////////////////// NODE ///////////////////////

/** @brief Builds a Node without connecting it to a tree.
 * @details This constructor allows the user to specify a physical range which
 * the node should be covering, and sets all references to parent or children
 * nodes to NULL. This constructor is useful to initialize the root of a tree
 * for example. Note another constructor is provided, which explicitely sets up
 * reference to a parent node, but the latter is private and will be called by
 * "wrapper" routines, and not the user directly.
 * @param[in] left,right,top,bottom Geometrical limits of the Node's domain.
 */
Node::Node( double left, double right, double top, double bottom )
    :Node(NULL,NULL)
{
    // Level, index, left and bottom are already initialized correctly.
    this->top_ = top; this->right_ = right;
    this->bottom_ = bottom; this->left_ = left;
}

/** @brief Constructor setting up the new Node's references in a tree
 * automatically.
 * @details This constructor accepts a reference to a parent node, as well as a
 * pointer to a particle that will be assigned to its center of mass attribute.
 * Multiple properties of the node, such as thee geometrical boundaries, the
 * level and index will be computed from the parent's data.
 * @param[in] parent Pointer to a parent node.
 * @param[in] part Pointer to a particle, which we assign to the center of mass.
 */
Node::Node( Node* parent, Particle* part )
    :parent_(parent), children_{NULL,NULL,NULL,NULL}, level_(0), index_(0), com_(part), left_(0.0), right_(0.0), top_(0.0), bottom_(0.0)
{
    // If parent exists, we can infer some properties (such as index, level,
    // boundaries) to the current node.
    if(parent_!=NULL) {
        level_ = parent_->level_+1;
        // Find which child we are and compute properties
        unsigned ic( parent_->GetQuadrant( part ) );
        // Shift the node's index by two bits compared to the parent,
        // then, add local index
        index_ = (parent_->index_ << 2) + ic;
        // Check child location relative to parent. Using the Z-curve
        // structure, if the zeroth bit of ic is one (ic=0,2), the node
        // is north. If the next bit of ic is one (ic=2,3), the node is
        // east.
        bool is_north( ic & 1U );
        bool is_east( (ic>>1) & 1U );
        left_   = parent_->left_   + 0.5*(double)(is_east)  *(parent_->right_-parent_->left_);
        right_  = parent_->right_  - 0.5*(double)(!is_east) *(parent_->right_-parent_->left_);
        bottom_ = parent_->bottom_ + 0.5*(double)(is_north) *(parent_->top_-parent_->bottom_);
        top_    = parent_->top_    - 0.5*(double)(!is_north)*(parent_->top_-parent_->bottom_);
    }
}

/** @brief Manages the destruction of Node objects and frees memory.
 * @details When a Node is deleted, all its initialized children are first
 * deleted. The memory is expected to have been allocated dynamically. Then, if
 * a Particle object is set for the Node, its identifier is compared to -1. If
 * the Particle is fictious (i.e. it represents the center of mass of multiple
 * particles in lower levels), it must have been allocated dynamically, and it
 * is deleted. Else, the Particle corresponds to an actual body that is being
 * simulated, and its memory was initialized by the IOManager class, and will be
 * destructed by the SConfig object it lives in.
 */
Node::~Node() {
    for( unsigned ic(0); ic<4; ic++) {
        if( this->children_[ic] != NULL ) {
            delete this->children_[ic];
            this->children_[ic] = NULL;
        }
        if( com_ != NULL ) {
            if( com_->id == -1 ){
                delete com_;
                com_ = NULL;
            }
        }
    }
}

/** Returns the parent node of the current object.
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

/** Answers true if the node is at the root of a tree (level zero).
 */
bool Node::IsRoot() const {
    return level_==0;
}

/** Answers true if all children of the node are null pointers.
 */
bool Node::IsLeaf() const {
    for( auto child : children_ )
        if( child != NULL )
            return false;
    return true;
}

/** Returns a pointer to the particle object of the node, that is, to its center
 * of mass.
 */
Particle* Node::GetParticle() const {
    return com_;
}

/** @brief Returns the center of mass of the node, which is gotten from the
 * local Particle object.
 */
std::array<double,2> Node::GetCenterOfMass() const {
    return this->com_->pos;
}

/** @brief Returns an array containing the node's geometrical bounds.
 */
std::array<double,4> Node::GetBounds() const {
    return std::array<double,4>{left_,right_,top_,bottom_};
}

/** Accessor to the current node's mass.
 * @returns The node's total mass.
 */
double Node::GetMass() const { return this->com_->mass; }

/** Returns the level of the current Node. The RootNode (likely a parent of the
 * current one) has level zero.
 */
long long Node::GetLevel() const { return level_; }

/** Returns the current Node's index in its local sub-tree.
 * @returns The Node's index.
 */
long long Node::GetIndex() const { return index_; }

/** Returns the width of the current node, that is, the length of one of its
 * edges (in this program, all nodes are squares).
 */
double Node::GetWidth() const {
    return right_ - left_;
}

/** @brief Returns true if the argmuent Particle is in the domain managed by the
 * caller node.
 * @param[in] p Pointer to a Particle object.
 * @returns Boolean indicating the Particle's belonging to the domain.
 */
bool Node::BelongsTo( Particle *p ) const {
    double xp(p->pos[0]); double yp(p->pos[1]);
    if( xp > this->left_ && xp < this->right_ &&
            yp > this->bottom_ && yp < this->top_ )
        return true;
    return false;
}

/** Returns true if no Particle is assigned to the current node.
 */
bool Node::IsEmpty() const {
    return (this->com_ == NULL );
}

/** For an input Particle, computes which quadrant of the current node (the
 * child index) would store this particle. Such a function is useful when
 * building a tree.
 */
unsigned Node::GetQuadrant( Particle* p ) const {
    unsigned idx(0);
    // If particle in top quadrants, toggle zeroth bit.
    if( p->pos[1] > bottom_ + 0.5*(top_-bottom_) )
        idx ^= 1U << 0;
    // If particle in right quadrants, toggle first bit.
    if( p->pos[0] > left_ + 0.5*(right_-left_) )
        idx ^= 1U << 1;
    return idx;
}

/** For a pointer to a node, computes the quadrant which would host it in the
 * caller node. The function operating on a particle pointer is used by defining
 * a Particle centered on the argument node's geometrical center, and passing it
 * to the other function.
 * @param[in] n Pointer to a node, supposedly located within the current node's
 * geometrical boundaries.
 * @returns Quadrant index with respect to the current node.
 */
unsigned Node::GetQuadrant( Node* n ) const {
    Particle tmp;
    tmp.pos = {0.5*(n->left_+n->right_), 0.5*(n->bottom_+n->top_)};
    return GetQuadrant( &tmp );
}

/** Returns the force applied along dim.
 * @param[in] dim Index of the dimension along which the force is requested
 * (0=x, 1=y).
 */
double Node::GetForce( unsigned dim ) const {
    /// @todo Delete this fucking mess.
    return com_->frc[dim];
}

/** Returns the number of initialized children of the caller.
 */
short Node::ChildrenCount() const {
    short nc(0);
    for( unsigned ic(0); ic<4; ic++ )
        if( this->children_[ic] != NULL )
            nc++;
    return nc;
}

/** Resets the center of mass of the current Node to default values.
 * This method is used during tree updates, before passing leaf information up
 * the tree.
 */
SError Node::ResetCenterOfMass() const {
    com_->mass = 0.0;
    com_->pos[0] = 0.0;
    com_->pos[1] = 0.0;
    com_->vel[0] = 0.0;
    com_->vel[1] = 0.0;
    ResetForces();
    return E_SUCCESS;
}

/** Resets the forces applied to a Node's center of mass.
 * This method is called by the Simulation::ComputeForce routine.
 */
SError Node::ResetForces() const {
    com_->frc[0] = 0.0;
    com_->frc[1] = 0.0;
    return E_SUCCESS;
}

/** Adds the gravitational force cause by the center of mass of a distant node
 * to the current node's center of mass.
 * The distant node's forces are not updated symmetrically.
 * @param[in] other Reference to a distant node.
 */
SError Node::Interact( const Node& other ) const {
    std::array<double,2> afrc = pp_force( *(this->com_), *(other.com_) );
    this->com_->frc[0] += afrc[0];
    this->com_->frc[1] += afrc[1];
    return E_SUCCESS;
}

/** @brief Prints Node information to an output stream.
 * @details This methods prints basic node information to the terminal,
 * including the node's memory address, its level, index, geometrical
 * boundaries, as well as the memory addresses of its children.
 * @note This overloading is thought for debugging purposes, and does not follow
 * strict conventions as to how the data is output. See IOManager for exporting
 * simulation results.
 * @param[in,out] os Output stream to which the information is printed.
 * @param[in] node Node which's information are to be printed.
 * @returns The input stream after instertion of the text to print.
 */
std::ostream& operator<<( std::ostream& os, const Node& node ) {
    os  << "[NODE " << &node << "] LVL=" << node.level_
        // << "/IDX=" << std::bitset<8*sizeof(long long)>(node.index_)
        << "/IDX=" << node.index_
        << " (L,R,B,T)=(" << node.left_ << "," << node.right_ << ","
        << node.bottom_ << "," << node.top_
        << ") CHILD=( ";
    for( unsigned ic(0); ic<4; ic++ )
        os << node.children_[ic] << " ";
    os << ")";
    return os;
}

/** Returns the squared distance between two nodes. This distance is computed
 * according to the node's center of mass, and not their geometrical center.
 */
double distance2( const Node& n1, const Node& n2 ) {
    const std::array<double,2>& pos1( n1.GetCenterOfMass() );
    const std::array<double,2>& pos2( n2.GetCenterOfMass() );
    return ( pos1[0]-pos2[0] )*( pos1[0]-pos2[0] ) + ( pos1[1]-pos2[1] )*( pos1[1]-pos2[1] );
}

/** Returns the distance between two nodes, by computing the square root of the
 * output of the function distance2. It is provided separately to avoid
 * performance issue.
 */
double distance( const Node& n1, const Node& n2 ) {
    return std::sqrt( distance2(n1,n2) );
}

/* ------------------------------------------- */

/** Builds a QuadTree object from a simulation configuration. It initializes the
 * root Node with the geometrical boundaries provided in the configuratio
 * object.
 * @param[in] config Reference to a configuration object containing the
 * simulation's parameters.
 */
QuadTree::QuadTree( const SConfig& config )
    :root_(new Node(0.0,config.dsize,config.dsize,0.0))
{}

/** Since the root node was allocated dynamically, destroy it if still set.
 * We note the destructors of the children of the root should be called by
 * Node::~Node.
 */
QuadTree::~QuadTree() {
    if( root_ != NULL ) {
        delete root_;
        root_ = NULL;
    }
}

/** @brief Tree browsing routine moving to the "right" of the tree.
 * @details Although the QuadTree class stores a pointer to the root of the
 * tree, it has little idea of the subsequent structure of this one. We provide
 * routines browsing the tree with the current function.
 * This methods looks for a Node with a higher index than the passed pointer,
 * `ptr`. The Node is first search among the parent of `ptr`, and if it cannot
 * be found (if `ptr` already has the maximal index in the branch), the routine
 * moves on level up in the tree and retries.
 * If no node can be found, which means we started with the highest node index
 * at any possible level, the method reaches the root node and returns NULL. If
 * a Node is returned, we note that its level might be different (higher) than
 * the one of `ptr`.
 * @param[in] ptr Pointer to the starting point of the movement routine.
 * @returns A pointer to the Node having a higher index than `ptr`, possibly at
 * a different level, or NULL if no Node can be found.
 */
Node* QuadTree::GetNext( Node* ptr ) const {
    // If we're at the Root, we want to return immediately.
    if( ptr->IsRoot() )
        return NULL;
    // Node is not last in quad tree decomposition. Parent has a child with
    // higher index (which we return, most common case).
    unsigned int idx(ptr->index_%4);
    for( unsigned uidx(idx+1); uidx<4; uidx++ ) {
        if( ptr->GetParent()->GetChild(uidx) != NULL ) {
            ptr = ptr->GetParent()->GetChild(uidx);
            return ptr;
        }
    }
    // Else, we need to search the next Node among the parents, get the current
    // node's.
    Node* nxt( ptr->GetParent() );
    // We stepped a level higher. Make sure we return the next node on the same
    // level.
    while( true ) {
        // If reached RootNode, return NULL; we come from the last child of this
        // RootNode and cannot look further.
        if( nxt->IsRoot() )
            return NULL;
        // This node's parent contains more interesting children, select the
        // next one (horizontal move in the tree).
        idx = nxt->index_%4;
        for( unsigned uidx(idx+1); uidx<4; uidx++ ) {
            if( nxt->GetParent()->GetChild(uidx) != NULL ) {
                nxt = nxt->GetParent()->GetChild(uidx);
                return nxt;
            }
        }
        // If this statement is reached, this means the investigated parent is
        // also the last child of the "upper" level, we need to search higher up
        // the tree.
        nxt = nxt->GetParent();
    }
}

/** Moves the passed pointer to the child node with lowest index. If no child is
 * found, the same pointer is returned.
 * @param[in] ptr Pointer to a node which we want to browse deeper.
 * @returns A pointer to the first available child of `ptr`.
 */
Node* QuadTree::GetDown( Node* ptr ) const {
    // Follow first child available
    for( unsigned ic(0); ic<4; ic++ ) {
        if( ptr->children_[ic] != NULL ) {
            ptr = ptr->children_[ic];
            break;
        }
    }
    return ptr;
}

/** @brief Looks for the next leaf of the tree from a given starting point.
 * @details From a given starting point, searches the closest leaf with greater
 * index. If the starting point, `ptr` is not a leaf, the method finds one by
 * recursively browsing lower levels. If the starting point is a leaf, the
 * method browses up the tree, until a parent with a child with a higher index
 * than the one of `ptr` is found. If no leaf can be found (we started from the
 * last one for example), the method reaches the root node and returns NULL.
 * @param[in] ptr Starting point of the tree search.
 * @returns A pointer to the next available leaf.
 */
Node* QuadTree::GetNextLeaf( Node* ptr ) const {
    if(ptr==NULL){
        return NULL;
    }
    Node* nxt(ptr);
    do{
        // Goes as deep as it can.
        if( not nxt->IsLeaf() ) {
            for( unsigned ic(0); ic<4; ic++ ) {
                if( nxt->GetChild(ic) != NULL ) {
                    nxt = nxt->GetChild(ic);
                    break;
                }
            }
        } else {
            // We reached a leaf. Need to go up until we can move again.
            bool go_up(true);
            while(go_up){
                // If we're at the root, cannot climb
                if(nxt->IsRoot()){
                    return NULL;
                }
                // Else get the parent.
                // Figure out our index from parent
                unsigned icur(nxt->GetParent()->GetQuadrant(nxt));
                nxt = nxt->GetParent();
                for( unsigned ic(icur+1); ic<4; ic++ ) {
                    if( nxt->GetChild(ic) != NULL ){
                        // Found a child from parent
                        nxt = nxt->GetChild(ic);
                        go_up = false;
                        break;
                    }
                }
                // Here, we're about to redo the while.
                // If go_up = false, we found a child to go to.
                // Else, the routine restarts and goes up one other level.
            }
        }
        // We only exit the routine if the candidate node is a leaf.
    }while(not nxt->IsLeaf());
    return nxt;
}

/** @brief Adds a Particle to the tree.
 * @details This method will look for an available child to insert the particle
 * in, as the node's center of mass. If a region where no Node is initialized is
 * found, the memory for this node is allocated dynamically, and the particle
 * assigned. If an occupied Node is found, the routine recursively splits this
 * Node (deepens the tree) until the particle we wish to insert, and the one
 * existing already fall into separate regions.
 * Each time a level is entered, the parent's center of mass is updated to
 * account for the particle that we will insert, such that upstream information
 * is already correct.
 * @param[in] p A pointer to the particle we wish to insert.
 * @returns An error code.
 * @todo This routine is not safe if a method tries to insert a copy of an
 * existing particle (it will deepen the tree forever).
 */
SError QuadTree::AddParticle( Particle* p ) const {
    Node* ptr(this->root_);
    while(true) {
        unsigned ic( ptr->GetQuadrant( p ) );
        if( ptr->IsLeaf() ) {
            // First off, check we're not at the Root (could mean tree is
            // empty!)
            if( ptr->com_ == NULL ) {
                ptr->com_ = p;
                break;
            }
            // We can insert right away, but we need also need to shift the
            // previous node down at the leaf level (else we lose a particle).
            // Problem:
            // The current particle, which needs to be shift down,
            // could end in the same level as the particle we're trying to
            // assign!
            // Solution:
            // Move the current node down, and let one more iteration run, so
            // another if statement is executed. The current node is replaced by
            // a fictious one.
            unsigned in( ptr->GetQuadrant( ptr->com_ ) );
            ptr->children_[in] = new Node( ptr, ptr->com_ );
            ptr->com_ = new Particle(*(ptr->com_));
        }
        if( ptr->children_[ic] == NULL ) {
            // The current node is not a leaf, but the quadrant we want to
            // insert the particle in is free, so we do it right away.
            // Increment the parent node with the mass of its new child.
            *(ptr->com_) += *p;
            ptr->children_[ic] = new Node( ptr, p );
            break;
        } else {
            // The current node is not a leaf, and the quadrant we need to
            // insert the particle in is already occupied. Move the pointer down
            // one level.
            // Increment the parent node with the mass of its new child.
            *(ptr->com_) += *p;
            ptr = ptr->children_[ic];
        }
    }
    return E_SUCCESS;
}

/** This method insert multiple particles, one after the other in the tree.
 */
SError QuadTree::AddParticle( std::vector<Particle*> p ) const {
    for( auto ap : p )
        this->AddParticle( ap );
    return E_SUCCESS;
}

/** Removes a particle from the tree. If a fictive particle is removed, its
 * memory is freed. Upstream information is updated for the parent nodes' center
 * of mass.
 * @param[in] ptr Pointer to the node containing the center of mass to delete
 * @returns An error code.
 * @note This method does NOT delete the Node, but leaves it empty. This choice
 * is motivated by the fact that deleting the Node would alter the tree
 * structure (sometimes massively) while the tree is browsed by other methods
 * (RemoveParticle is called by Simulation::UpdateTree). We leave the deletion
 * of unecessary nodes to a more robust routine, QuadTree::PruneNode.
 */
SError QuadTree::RemoveParticle( Node* ptr ) const {
    Node *node( ptr->GetParent() );
    // Remove the particle's upstream influence in higher levels.
    while( node != NULL ){
        *(node->com_) -= *(ptr->com_);
        node = node->GetParent();
    }
    // If the Particle is a fictive one, delete it.
    /// @todo This is definitely not a good idea.
    if( ptr->com_->id == -1 ) {
        delete ptr->com_;
    }
    ptr->com_ = NULL;
    return E_SUCCESS;
}

/** @brief Deletes a Node from the tree.
 * @details Frees the memory associated with the argument Node pointer. The
 * parent node is then checked, if it contains only one child, this means the
 * structure of the tree can be simplified (levels can be removed). The
 * remaining child is found, and "dragged" to higher levels if it is itself a
 * leaf (if not, the deeper structure is still needed).
 * @param[in] ptr A pointer to the Node we wish to delete.
 * @returns An error code.
 */
SError QuadTree::PruneNode( Node* ptr ) const {
    // Prune the Node like demanded
    Node* n( ptr->GetParent() );
    unsigned ic( ptr->GetIndex() % 4 );
    delete n->children_[ic];
    n->children_[ic] = NULL;
    while( n != NULL ) {
        if( n->ChildrenCount() != 1 ) {
            break;
        } else {
            // Only one children remains in the parent of the pruned node. This
            // child can be merged with its parent.
            // One needs to care about the particles within the Nodes being
            // fictive or not (id==-1) when performing this operation, so that
            // that the output remains coherent.
            Node* child( NULL );
            // Find the remaining child.
            for( ic=0; ic<4; ic++ ) {
                if( n->children_[ic] != NULL ) {
                    child = n->children_[ic];
                    break;
                }
            }
            // Check if the remaining node is  a leaf. If it is, we may drag it
            // up one level. If not, multiple particles are enclosed close to
            // one another at depper levels, so do not touch anything.
            if( child->IsLeaf() ) {
                // Swap contained particle pointers (move the leaf particle to the
                // parent).
                std::swap( n->com_, child->com_ );
                // Delete the remaining only child
                delete n->children_[ic];
                n->children_[ic] = NULL;
                // Get the parent and perform the same checks again
                n = n->GetParent();
            } else {
                break;
            }
        }
    }
    return E_SUCCESS;
}

/** @brief Overloaded operator for debugging purposes.
 * @details This method includes the most
 * complete tree browsing routine available in this class, QuadTree::GetNext and
 * QuadTree::GetNextLeaf are just shortened versions of the one available here.
 * @param[in,out] os Output stream to which data is written.
 * @param[in] tree Reference to the QuadTree object to print.
 * @returns An output stream to which information has been written.
 * @note This method is unsuitable for post-analysis-ready output. Please refer
 * to IOManager for this purpose.
 */
std::ostream& operator<<( std::ostream& os, const QuadTree& tree ) {
    os << "[TREE " << &tree << "] root_=" << tree.root_ << std::endl;
    Node* ptr(tree.root_);
    os << "    " << *ptr << std::endl
       << "          " << *(ptr->com_) << std::endl;
    while(true) {
        // Goes as deep as it can, prints all nodes when stepping down.
        if( not ptr->IsLeaf() ) {
            for( unsigned ic(0); ic<4; ic++ ) {
                if( ptr->GetChild(ic) != NULL ) {
                    ptr = ptr->GetChild(ic);
                    os << "    " << *ptr << std::endl
                       << "          " << *(ptr->com_) << std::endl;
                    break;
                }
            }
        } else {
            // We reached a leaf. Need to go up until we can move again.
            bool go_up(true);
            while(go_up){
                // If we're at the root, cannot climb
                if(ptr->IsRoot()){
                    return os;
                }
                // Else get the parent.
                // Figure out our index from parent
                unsigned icur(ptr->GetParent()->GetQuadrant(ptr));
                ptr = ptr->GetParent();
                for( unsigned ic(icur+1); ic<4; ic++ ) {
                    if( ptr->GetChild(ic) != NULL ){
                        // Found a child from parent
                        ptr = ptr->GetChild(ic);
                        os << "    " << *ptr << std::endl
                           << "          " << *(ptr->com_) << std::endl;
                        go_up = false;
                        break;
                    }
                }
                // Here, we're about to redo the while.
                // If go_up = false, we found a child to go to.
                // Else, the routine restarts and goes up one other level.
            }
        }
    }
    return os;
}
