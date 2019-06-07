#ifndef SURMISE_SROUTINES_HPP_
#define SURMISE_SROUTINES_HPP_

#include <mpi.h>
#include "stypes.hpp"
#include "timer.hpp"
#include "iomanager.hpp"
#include "nodes.hpp"
#include "serrors.hpp"

extern MPI_Datatype MPI_Particle;

/** Wrapper around the datatypes defined in this simulation. It runs a N-body
 * gravity simulation, manages the outputs, and stuff like that.
 */
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
        /**
         */
        Timer* timer_;
        /** Constant reference to the simulation configuration produced by the
         * IOManager class. It contains all particles of the system, dynamically
         * allocated by the input routines.*/
        const SConfig& conf_;
        /** Pointer to a class containing a quad-tree. This tree will be by the
         * Barnes-Hut algorithm.*/
        QuadTree* tree_;
        /** List of boolean flags determining if the current simulation instance
         * is responsible of an existing particle or not.
         */
        bool* update_list_;
};


#endif // SURMISE_SROUTINES_HPP_
