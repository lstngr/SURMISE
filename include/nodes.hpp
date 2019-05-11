/** @file nodes.hpp
 * @brief Defines different types of nodes (root, normal, leaf) that will be
 * used for domain decomposition.
 */
#ifndef SURMISE_NODES_HPP_
#define SURMISE_NODES_HPP_

#include "stypes.hpp"

/** @class Node
 * @brief Generic Node in the decomposition tree.
 * @details It contains information about its parent domain, children, and has
 * relevant methods performing FMM routines, particle distribution routines and
 * so on.
 */
class Node {
    public:
        Node(Node *parent = NULL );
        virtual ~Node(){};
        Node* GetParent() const;
        Node* GetChild(short int child_idx) const;
        virtual int Decompose( Particle* parts );
        virtual int TimeEvolution( double dt );
        unsigned int GetLevel() const;
    protected:
    private:
        /** Pointer to parent Node*/
        Node *parent_;
        /**Array of pointers to children node. As we work in the two-dimensional
         space, the size is preallocated to a quad-tree structure.*/
        Node *children_[4];
        /** Level of the Node. The RootNode has level 0.*/
        unsigned int level_;
        ///@{
        /// Geometrical limits of the node's coverage (of physical space)
        float xstart_, xend_, ystart_, yend_;
        ///@}
};

/** @brief Node at the root of the tree decomposition of the system.
 * @details This object is instanciated at the start of the simulation and
 * contains various methods to manage, run and export the simulation. It
 * represents the interface called from the main function.
 */
class RootNode : public Node {
    public:
        RootNode(SConfig& conf);
        int Run();
        int Decompose( Particle* parts ) override;
        int TimeEvolution( double dt ) override;
    protected:
    private:
        /** Configuration of the simulation. We note it contains the relevant
         * algorithmic parameters, but also the Particle's to be simulated.*/
        SConfig conf_;
};

/** @class LeafNode
 * @brief Node located at the bottom of the tree.
 * @details It handles the concrete evolution and computation of the system,
 * namely since it handles the Particle objects directly (direct force
 * summation, lowest level FMM, time evolution).
 */
class LeafNode : public Node {
    public:
        LeafNode( Node *parent, Particle* parts );
        int Decompose( Particle* parts ) override;
        int TimeEvolution( double dt ) override;
    protected:
    private:
        /** Stores the Particle located in the node. They can be used for direct
         * or nearest-neighbours computation.*/
        Particle *subparts_;
};

#endif //SURMISE_NODES_HPP_
