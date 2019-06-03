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
    /* ---------------------------- */
    // Define the MPI Particle type
    /* ---------------------------- */
    // Inspired by the example 4.17 of the MPI3.1 manual
    MPI_Datatype MPI_Particle, MPI_Particle_NR;
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
    MPI_Type_create_struct( 2, blocklen, disp, ptypes, &MPI_Particle_NR );
    // NOTE: Doubles and longs have different alignement (double have 8 bytes,
    // while long only four), such that the compiler might add a four byte
    // padding at the end of the struct. We resize it (will also allow to send
    // multiple particles at once). 
    int sz;
    MPI_Aint lb,ext;
    MPI_Type_size( MPI_Particle_NR, &sz );
    MPI_Type_get_extent( MPI_Particle_NR, &lb, &ext );
    if(sz != sizeof( Particle )) {
        MPI_Type_create_resized( MPI_Particle_NR,
                0,
                (MPI_Aint)sizeof(Particle),
                &MPI_Particle);
    } else {
        MPI_Particle = MPI_Particle_NR;
    }
    MPI_Type_commit( &MPI_Particle );
    /* ---------------------------- */
    // Run simulation
    IOManager io( argv[1] );
    Simulation sim( io );
    sim.Run();
    /* ---------------------------- */
    // Free MPI_Particle type.
    MPI_Type_free( &MPI_Particle );
    MPI_Finalize();
    return 0;
}
