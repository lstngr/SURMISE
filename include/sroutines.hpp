/**
 * @file Something
 */

#ifndef SURMISE_SROUTINES_HPP_
#define SURMISE_SROUTINES_HPP_

#include "stypes.hpp"
#include "iomanager.hpp"
#include "nodes.hpp"
#include "serrors.hpp"

class Simulation {
    public:
        Simulation( IOManager& io );
        ~Simulation();
        SError Run();
    private:
        friend class IOManager;
        SError BuildTree();
        SError UpdateTree() const;
        SError ComputeForces() const;
        SError TimeEvolution() const;
        IOManager& io_;
        const SConfig& conf_;
        QuadTree* tree_;
};


#endif // SURMISE_SROUTINES_HPP_
