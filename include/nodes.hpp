/** @file nodes.hpp
 * @brief Defines different types of nodes (root, normal, leaf) that will be
 * used for domain decomposition.
 */
#ifndef SURMISE_NODES_HPP_
#define SURMISE_NODES_HPP_

#include "stypes.hpp"

/** @brief Generic Node in the decomposition tree.
 * @details It contains information about its parent domain, children, and has
 * relevant methods performing FMM routines, particle distribution routines and
 * so on.
 * @var parent_ @private Pointer to parent Node
 * @var children_ @private Array of pointers to children node. As we work
 * in the two-dimensional space, the size is preallocated to a quad-tree
 * structure.
 * @var level_ @private Level of the Node. The RootNode has level 0.
 */
class Node {
    public:
        Node(Node *parent = NULL );
        virtual ~Node(){};
        Node* GetParent() const;
        Node* GetChild(short int child_idx) const;
        int Decompose();
        unsigned int GetLevel() const;
    protected:
    private:
        Node *parent_;
        Node *children_[4];
        unsigned int level_;
};

class RootNode : public Node {
    public:
        RootNode(SConfig& conf);
        int Run();
    protected:
    private:
        SConfig conf_;
};

class LeafNode : public Node {
    public:
    protected:
    private:
        Particle *subparts_;
};

#endif //SURMISE_NODES_HPP_
