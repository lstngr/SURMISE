/** @file nodes.hpp
 * @brief Defines different types of nodes (root, normal, leaf) that will be
 * used for domain decomposition.
 */
#ifndef SURMISE_NODES_HPP_
#define SURMISE_NODES_HPP_

#include <array>
#include <vector>
#include "stypes.hpp"
#include "serrors.hpp"

/** Enumeration associating a direction and an int. These directions are to be
 * interpreted in the 2D plane of the simulation, and are used by the methods
 * lookng for nearest neighbours.
 */
enum ZDIR {
    UP=1, DOWN=-1, LEFT=-2, RIGHT=2
};

/** @class Node
 * @brief Generic Node in the decomposition tree.
 * @details It contains information about its parent domain, children, and has
 * relevant methods performing FMM routines, particle distribution routines and
 * so on.
 */
class Node {
    public:
        Node(Node *parent = NULL, unsigned int idx = 0 );
        virtual ~Node();
        Node* GetParent() const;
        Node* GetChild(short int child_idx) const;
        Node* GetNext() const;
        Node* Move( ZDIR direction ) const;
        bool IsNN(const Node* other) const;
        std::array<Node*,8>  GetNN() const;
        std::array<Node*,27> GetIN() const;
        SError InitChild( int child );
        SError InitLeaf( int leaf );
        SError InitAllChildren();
        SError InitAllLeafs();
        std::array<double,2> Center() const;
        std::array<double,4> GetBounds() const;
        SError SetBounds(double x1, double x2, double y1, double y2);
        virtual SError Decompose( std::vector<Particle*>& parts );
        virtual SError TimeEvolution( double dt );
        virtual SError TimeEvolutionMasses( double dt );
        unsigned int GetLevel() const;
        unsigned int GetIndex() const;
        bool BelongsTo( Particle *p ) const;
        virtual SError Reassign( Particle* p );
        virtual SError AddParticle( Particle* p);
        static unsigned int maximum_level; /*!< Maximum level of the tree
                                    @todo Protect scope, check MPI feasability.*/
        SError Interact( const Node* other );
        virtual SError GatherMasses();
        double GetMass() const;
        SError SetMass(double m);
        SError ResetForces();
        SError AddForce(const std::array<double,2>& upstream_force);
        std::array<double,2> GetForce() const;
        /** Returns an empty vector of particles. @todo Ugly AF.*/
        virtual std::vector<Particle*> GetParticles() const{return std::vector<Particle*>{};}
    protected:
    private:
        /** Pointer to parent Node*/
        Node *parent_;
        /**Array of pointers to children node. As we work in the two-dimensional
         space, the size is preallocated to a quad-tree structure.*/
        Node *children_[4];
        /** Level of the Node. The RootNode has level 0.*/
        unsigned int level_;
        /** Index of the Node within its parent's children subtree. */
        unsigned int idx_;
        /// Geometrical limits of the node's coverage (of physical space)
        double xybnds[4];
        /** Total mass of the Node.*/
        double mass;
        /**Force applied on the Node's _geometrical_ center.*/
        std::array<double,2> force;
};

/** @brief Node at the root of the tree decomposition of the system.
 * @details This object is instanciated at the start of the simulation and
 * contains various methods to manage, run and export the simulation. It
 * represents the interface called from the main function.
 */
class RootNode : public Node {
    public:
        RootNode(SConfig& conf);
        SError Run();
        SError Decompose( std::vector<Particle*>& parts ) override;
        SError TimeEvolution( double dt ) override;
        SError Reassign( Particle* p ) override;
    protected:
    private:
        /** Configuration of the simulation. We note it contains the relevant
         * algorithmic parameters, but also the Particle's to be simulated.*/
        SConfig &conf_;
};

/** @class LeafNode
 * @brief Node located at the bottom of the tree.
 * @details It handles the concrete evolution and computation of the system,
 * namely since it handles the Particle objects directly (direct force
 * summation, lowest level FMM, time evolution).
 */
class LeafNode : public Node {
    public:
        LeafNode( Node *parent, unsigned int child );
        LeafNode( Node *parent, unsigned int child, std::vector<Particle*> &parts );
        SError Decompose( std::vector<Particle*>& parts ) override;
        SError TimeEvolution( double dt ) override;
        SError TimeEvolutionMasses( double dt ) override;
        SError Reassign( Particle* p ) override;
        SError AddParticle( Particle* p) override;
        SError GatherMasses() override;
        std::vector<Particle*> GetParticles() const override;
    protected:
    private:
        /** Stores the Particle located in the node. They can be used for direct
         * or nearest-neighbours computation.*/
        std::vector<Particle*> subparts_;
};

#endif //SURMISE_NODES_HPP_
