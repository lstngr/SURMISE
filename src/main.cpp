/** \file main.cpp
 * @brief This file launches the simulation.
 * @details In fact this is far more complicated than it looks like!
 */
#include <iostream>
#include <mpi.h>
#include "stypes.hpp"
#include "sinput.hpp"
#include "nodes.hpp"
#include "sroutines.hpp"

int main( int argc, char** argv ){
    MPI_Init( &argc, &argv );
    make_mpi_types();
    /* ---------------------------- */
    // Run simulation
    IOManager io( argv[1] );
    Simulation sim( io );
    MPI_Barrier( MPI_COMM_WORLD );
    sim.Run();
    /* ---------------------------- */
    // Free MPI_Particle type.
    free_mpi_types();
    MPI_Finalize();
    return 0;
}
