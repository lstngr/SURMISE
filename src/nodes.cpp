#include <iostream>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "nodes.hpp"

//////////////////// NODE ///////////////////////

Node::Node( double left, double right, double top, double bottom )
    :Node(NULL,NULL)
{
    // Level, index, left and bottom are already initialized correctly.
    this->top_ = top; this->right_ = right;
    this->bottom_ = bottom; this->left_ = left;
}

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
        // We found the node, break iteration
    }
}

/** @brief Manages the destruction of Node objects and memory management.
 * @details Moooar
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

bool Node::IsRoot() const {
    return level_==0;
}

bool Node::IsLeaf() const {
    for( auto child : children_ )
        if( child != NULL )
            return false;
    return true;
}

Particle* Node::GetParticle() const {
    return com_;
}

/** @brief Computes the geometrical center of the calling node from its
 * boundaries.
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

double Node::GetWidth() const {
    return right_ - left_;
}

/** @brief Returns true if the argmuent Particle is in the domain
 * @details Many more things
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

unsigned Node::GetQuadrant( Node* n ) const {
    Particle tmp;
    tmp.pos = {0.5*(n->left_+n->right_), 0.5*(n->bottom_+n->top_)};
    return GetQuadrant( &tmp );
}

double Node::GetForce( unsigned dim ) const {
    /// @todo Delete this fucking mess.
    return com_->frc[dim];
}

SError Node::ResetForces() const {
    com_->frc[0] = 0.0;
    com_->frc[1] = 0.0;
    return E_SUCCESS;
}

SError Node::Interact( const Node& other ) const {
    std::array<double,2> afrc = pp_force( *(this->com_), *(other.com_) );
    this->com_->frc[0] += afrc[0];
    this->com_->frc[1] += afrc[1];
    return E_SUCCESS;
}

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

double distance( const Node& n1, const Node& n2 ) {
    const std::array<double,2>& pos1( n1.GetCenterOfMass() );
    const std::array<double,2>& pos2( n2.GetCenterOfMass() );
    return std::sqrt( ( pos1[0]-pos2[0] )*( pos1[0]-pos2[0] ) + ( pos1[1]-pos2[1] )*( pos1[1]-pos2[1] ) );
}

/* ------------------------------------------- */

QuadTree::QuadTree( const SConfig& config )
    :root_(new Node(0.0,config.dsize,config.dsize,0.0))
{}

QuadTree::~QuadTree() {
    if( root_ != NULL ) {
        delete root_;
        root_ = NULL;
    }
}

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

SError QuadTree::AddParticle( std::vector<Particle*> p ) const {
    for( auto ap : p )
        this->AddParticle( ap );
    return E_SUCCESS;
}

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
