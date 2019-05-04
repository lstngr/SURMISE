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
    return 0;
}
