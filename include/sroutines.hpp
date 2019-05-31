/**
 * @file Something
 */

#ifndef SURMISE_SROUTINES_HPP_
#define SURMISE_SROUTINES_HPP_

#include "stypes.hpp"
#include "nodes.hpp"
#include "serrors.hpp"

class Simulation {
    public:
        Simulation( const SConfig& conf ) :conf_(conf){}
        ~Simulation();
        SError Run();
    private:
        friend class IOManager;
        SError BuildTree();
        SError UpdateTree() const;
        SError ComputeForces() const;
        SError TimeEvolution() const;
        const SConfig& conf_;
        QuadTree* tree_;
};


#endif // SURMISE_SROUTINES_HPP_
