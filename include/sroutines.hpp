/**
 * @file Something
 */

#ifndef SURMISE_SROUTINES_HPP_
#define SURMISE_SROUTINES_HPP_

#include <mpi.h>
#include "stypes.hpp"
#include "iomanager.hpp"
#include "nodes.hpp"
#include "serrors.hpp"

extern MPI_Datatype MPI_Particle;

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
        /** Reference to an input/output manager for configuration reading and
         * data output.*/
        IOManager& io_;
        /** Constant reference to the simulation configuration produced by the
         * IOManager class. It contains all particles of the system, dynamically
         * allocated by the input routines.*/
        const SConfig& conf_;
        /** Pointer to a class containing a quad-tree. This tree will be by the
         * Barnes-Hut algorithm.*/
        QuadTree* tree_;
        bool* update_list_;
};


#endif // SURMISE_SROUTINES_HPP_
