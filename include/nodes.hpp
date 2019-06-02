/** @file nodes.hpp
 * @brief Defines different types of nodes (root, normal, leaf) that will be
 * used for domain decomposition.
 */
#ifndef SURMISE_NODES_HPP_
#define SURMISE_NODES_HPP_

#include <array>
#include <vector>
#include <iostream>
#include "stypes.hpp"
#include "serrors.hpp"

/** Enumeration associating a direction and an int. These directions are to be
 * interpreted in the 2D plane of the simulation, and are used by the methods
 * lookng for nearest neighbours.
 */
enum ZDIR {
    UP=1, DOWN=-1, LEFT=-2, RIGHT=2
};

class QuadTree;

/** @class Node
 * @brief Generic Node in the decomposition tree.
 * @details It contains information about its parent domain, children, and has
 * relevant methods performing FMM routines, particle distribution routines and
 * so on.
 */
class Node {
    public:
        Node( double left, double right, double top, double bottom );
        virtual ~Node();
        Node* GetParent() const;
        Node* GetChild(short int child_idx) const;
        bool IsRoot() const;
        bool IsLeaf() const;
        Particle* GetParticle() const;
        std::array<double,2> GetCenterOfMass() const;
        std::array<double,4> GetBounds() const;
        long long GetLevel() const;
        long long GetIndex() const;
        double GetWidth() const;
        double GetMass() const;
        bool BelongsTo( Particle *p ) const;
        bool IsEmpty() const;
        unsigned GetQuadrant( Particle* p ) const;
        unsigned GetQuadrant( Node* n ) const;
        double GetForce( unsigned dim ) const;
        short ChildrenCount() const;
        SError ResetCenterOfMass() const;
        SError ResetForces() const;
        SError Interact( const Node& other ) const;
        friend class QuadTree;
        friend class IOManager;
        friend std::ostream& operator<<( std::ostream& os, const Node& node );
        friend std::ostream& operator<<( std::ostream& os, const QuadTree& tree );
    private:
        Node(Node *parent, Particle *part);
        /** Pointer to parent Node*/
        Node *parent_;
        /**Array of pointers to children node. As we work in the two-dimensional
         space, the size is preallocated to a quad-tree structure.*/
        Node *children_[4];
        unsigned long long level_;
        unsigned long long index_;
        Particle* com_;
        /// Geometrical limits of the node's coverage (of physical space)
        double left_,right_,top_,bottom_;
};

double distance( const Node& n1, const Node& n2 );
double distance2( const Node& n1, const Node& n2 );

class QuadTree {
    public:
        QuadTree( const SConfig& config );
        ~QuadTree();
        Node* GetRoot() const { return root_; }
        Node* GetNext( Node* ptr ) const;
        Node* GetDown( Node* ptr ) const;
        Node* GetNextLeaf( Node* ptr ) const;
        SError AddParticle( Particle* p ) const;
        SError AddParticle( std::vector<Particle*> p ) const;
        SError RemoveParticle( Node* ptr ) const;
        SError PruneNode( Node* ptr ) const;
        friend std::ostream& operator<<( std::ostream& os, const QuadTree& tree );
    private:
        Node* root_;
};

#endif //SURMISE_NODES_HPP_
