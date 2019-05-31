/** \file main.cpp
 * @brief This file launches the simulation.
 * @details In fact this is far more complicated than it looks like!
 */
#include <iostream>
#include "sinput.hpp"
#include "nodes.hpp"
#include "sroutines.hpp"

int main( int argc, char** argv ){
    IOManager io( argv[1] );
    Simulation sim( io );
    sim.Run();
    return 0;
}
