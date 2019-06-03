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
    /* ---------------------------- */
    // Define the MPI Particle type
    /* ---------------------------- */
    // Inspired by the example 4.17 of the MPI3.1 manual
    MPI_Datatype MPI_Particle;
    // Define structure types and block lengths
    MPI_Datatype ptypes[5]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_LONG};
    int blocklen[5]={2,2,2,1,1};
    // Compute displacements
    MPI_Aint disp[5];
    disp[0] = offsetof( Particle, pos );
    disp[1] = offsetof( Particle, vel );
    disp[2] = offsetof( Particle, frc );
    disp[3] = offsetof( Particle, mass);
    disp[4] = offsetof( Particle, id  );
    MPI_Type_create_struct( 2, blocklen, disp, ptypes, &MPI_Particle );
    MPI_Type_commit( &MPI_Particle );
    /* ---------------------------- */
    // Run simulation
    IOManager io( argv[1] );
    Simulation sim( io );
    sim.Run();
    MPI_Finalize();
    return 0;
}
