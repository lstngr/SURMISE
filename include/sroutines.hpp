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
        Simulation( SConfig& conf ) :conf_(conf){}
        ~Simulation();
        SError Run();
    private:
        SError BuildTree();
        SError UpdateTree();
        SError ComputeForces();
        SError TimeEvolution();
        SConfig conf_;
        QuadTree* tree_;
};


#endif // SURMISE_SROUTINES_HPP_
