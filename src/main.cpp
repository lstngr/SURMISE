/** \file main.cpp
 * @brief This file launches the simulation.
 * @details In fact this is far more complicated than it looks like!
 */
#include <iostream>
#include "sinput.hpp"
//#include "sroutines.hpp"

int main( int argc, char** argv ){
    int errc(0);
    SConfig sim_conf(ReadInputs(argc,argv,errc));
    for(size_t ip(0); ip<sim_conf.npart; ip++){
        std::cout << sim_conf.parts[ip].pos[1] << std::endl;
    }
    delete[] sim_conf.parts; //TODO - Move to destructor
    return 0;
}
