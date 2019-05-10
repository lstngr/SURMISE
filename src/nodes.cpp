#include "nodes.hpp"

//////////////////// NODE ///////////////////////

/** @brief Creates a generic node.
 * @param[in] parent Pointer to a parent node. If none is provided, the pointer
 * is left null and the node is interpreted to be at the root.
 */
Node::Node( Node *parent )
    :parent_(parent)
{
    if( parent_ != NULL )
        this->level_ = parent_->GetLevel() + 1;
    else
        this->level_ = 0;
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

//////////////////// ROOTNODE ///////////////////////
RootNode::RootNode( SConfig& cf )
    :conf_(cf)
{}

//////////////////// LEAFNODE ///////////////////////

