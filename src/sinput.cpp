#include <string>
#include <iostream>
#include <fstream>
#include "sinput.hpp"
#include "config_file.hpp"

SConfig ReadInputs( int argc, char** argv, int& error_code ) {
    // Declare simulation configuration object
    if(argc != 3){
        throw;
    }
    SConfig sim_conf;
    /// @todo - Error management
    error_code = ReadInputFile(sim_conf, argv[1], argv[2]);
    return sim_conf;
}

int ReadInputFile( SConfig& sim_conf, const std::string& filename, const std::string& filename_init ) {
    ConfigFile cfg(filename);
    // Load basic parameters
    sim_conf.dsize      = cfg.get<double>("domsize");
    sim_conf.npart      = cfg.get<int>("npart");
    sim_conf.dt         = cfg.get<double>("tevol_dt");
    sim_conf.epsilon    = cfg.get<double>("epsilon");
    sim_conf.max_iter   = cfg.get<int>("max_iter");
    sim_conf.max_wtime  = cfg.get<double>("walltime");
    sim_conf.extra_time = cfg.get<double>("extratime");
    // Load particles
    // Particle array is dynamically allocated and should match input file size.
    sim_conf.parts.reserve(sim_conf.npart);
    std::ifstream ifsp( filename_init );
    for( size_t ip(0); ip < sim_conf.npart; ip++  ) {
        sim_conf.parts.push_back( new Particle );
        sim_conf.parts.back()->id = ip;
        ifsp >> sim_conf.parts.back()->pos[0]
             >> sim_conf.parts.back()->pos[1]
             >> sim_conf.parts.back()->vel[0]
             >> sim_conf.parts.back()->vel[1]
             >> sim_conf.parts.back()->mass;
        sim_conf.parts.back()->frc[0] = 0.0;
        sim_conf.parts.back()->frc[1] = 0.0;
    }
    ifsp.close();
    return 0;
}
