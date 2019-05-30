#include <iostream>
#include <iomanip>
#include <cmath>
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
Node* QuadTree::GetNext( Node* ptr ) const {
    // If we're at the Root, we want to return immediately.
    if( ptr->level_ == 0 )
        return NULL;
    // Node is not last in quad tree decomposition. Parent has a child with
    // higher index (which we return, most common case).
    unsigned int idx(ptr->index_%4);
    if( idx<3 )
        return ptr->GetParent()->GetChild( idx + 1 );
    // Else, we need to search the next Node among the parents, get the current
    // node's.
    Node* nxt( ptr->GetParent() );
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
            while( nxt->GetLevel() != ptr->level_ )
                nxt = nxt->GetChild( 0 );
            return nxt;
        }
        // If this statement is reached, this means the investigated parent is
        // also the last child of the "upper" level, we need to search higher up
        // the tree.
        nxt = nxt->GetParent();
    }
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

std::ostream& operator<<( std::ostream& os, const Node& node ) {
    os << "[NODE " << &node << "] LVL=" << node.level_ << "/IDX=" << node.index_
        << " (L,R,B,T)=(" << node.left_ << "," << node.right_ << ","
        << node.bottom_ << "," << node.top_
        << ") CHILD=( ";
    for( unsigned ic(0); ic<4; ic++ )
        os << node.children_[ic] << " ";
    os << ")";
    return os;
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

SError QuadTree::AddParticle( Particle* p ){
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
            // another if statement is executed.
            unsigned in( ptr->GetQuadrant( ptr->com_ ) );
            ptr->children_[in] = new Node( ptr, ptr->com_ );
        }
        if( ptr->children_[ic] == NULL ) {
            // The current node is not a leaf, but the quadrant we want to
            // insert the particle in is free, so we do it right away.
            ptr->children_[ic] = new Node( ptr, p );
            break;
        } else {
            // The current node is not a leaf, and the quadrant we need to
            // insert the particle in is already occupied. Move the pointer down
            // one level.
            ptr = ptr->children_[ic];
        }
    }
    return E_SUCCESS;
}

SError QuadTree::AddParticle( std::vector<Particle*> p ) {
    for( auto ap : p )
        this->AddParticle( ap );
    return E_SUCCESS;
}

std::ostream& operator<<( std::ostream& os, const QuadTree& tree ) {
    os << "[TREE " << &tree << "] root_=" << tree.root_;
    return os;
}
