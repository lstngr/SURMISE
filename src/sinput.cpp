#include <unistd.h>
#include <string>
#include "sinput.hpp"
#include "config_file.hpp"

SConfig ReadInputs( int argc, char** argv, int& error_code ) {
    // Declare simulation configuration object
    SConfig sim_conf;
    // Declared in unistd.h (Part if the getopt() routines)
    opterr = 0;
    // Argument counter
    int c;
    /** Scanned command line arguments:
     * - 'f': Input file containing variables to set at start
     */
    while ((c = getopt (argc, argv, "f:")) != -1) {
        switch (c) {
            case 'f':
                error_code = ReadInputFile(sim_conf,optarg);
                break;
            case '?':
                if (optopt == 'f')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint(optopt))
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                error_code=3;
            default:
                std::cout << "Coucou" << std::endl;
        }
    }
    return sim_conf;
}

int ReadInputFile( SConfig& sim_conf, const std::string& filename ) {
    ConfigFile cfg(filename);
    // Load basic parameters
    sim_conf.dsize      = cfg.get<double>("domsize");
    sim_conf.npart      = cfg.get<int>("npart");
    sim_conf.dt         = cfg.get<double>("tevol_dt");
    sim_conf.epsilon    = cfg.get<double>("epsilon");
    sim_conf.max_iter   = cfg.get<int>("max_iter");
    sim_conf.max_wtime  = cfg.get<double>("max_wtime");
    sim_conf.extra_time = cfg.get<double>("extra_time");
    // Load initial conditions
    std::string dxlocs(cfg.get<std::string>("pos_x_distrib"));
    std::string dylocs(cfg.get<std::string>("pos_x_distrib"));
    std::string dxvels(cfg.get<std::string>("vel_x_distrib"));
    std::string dyvels(cfg.get<std::string>("vel_x_distrib"));
    std::string dmass (cfg.get<std::string>("mass_distrib"));
    // Position profile
    if( dxlocs=="uniform" ){
        sim_conf.distr_xlocs = uniform;
        sim_conf.unifo_xmin = cfg.get<double>("x_uniform_min");
        sim_conf.unifo_xmax = cfg.get<double>("x_uniform_max");
    } else if( dxlocs=="gauss" ){
        sim_conf.distr_xlocs = gauss;
        sim_conf.gauss_xmean = cfg.get<double>("x_gauss_mean");
        sim_conf.gauss_xstdd = cfg.get<double>("x_gauss_stdd");
    } else {
        std::cerr << "Unknown initial x-position profile" << std::endl;
    }
    if( dylocs=="uniform" ){
        sim_conf.distr_ylocs = uniform;
        sim_conf.unifo_ymin = cfg.get<double>("y_uniform_min");
        sim_conf.unifo_ymax = cfg.get<double>("y_uniform_max");
    } else if( dylocs=="gauss" ){
        sim_conf.distr_ylocs = gauss;
        sim_conf.gauss_ymean = cfg.get<double>("y_gauss_mean");
        sim_conf.gauss_ystdd = cfg.get<double>("y_gauss_stdd");
    } else {
        std::cerr << "Unknown initial y-position profile" << std::endl;
    }
    // Velocity profile
    if( dxvels=="uniform" ){
        sim_conf.distr_xvels = uniform;
        sim_conf.unifo_vxmin = cfg.get<double>("vx_uniform_min");
        sim_conf.unifo_vxmax = cfg.get<double>("vx_uniform_max");
    } else if( dxvels=="gauss" ){
        sim_conf.distr_xvels = gauss;
        sim_conf.gauss_vxmean = cfg.get<double>("vx_gauss_mean");
        sim_conf.gauss_vxstdd = cfg.get<double>("vx_gauss_stdd");
    } else {
        std::cerr << "Unknown initial x-velocity profile" << std::endl;
    }
    if( dyvels=="uniform" ){
        sim_conf.distr_yvels = uniform;
        sim_conf.unifo_vymin = cfg.get<double>("vy_uniform_min");
        sim_conf.unifo_vymax = cfg.get<double>("vy_uniform_max");
    } else if( dyvels=="gauss" ){
        sim_conf.distr_yvels = gauss;
        sim_conf.gauss_vymean = cfg.get<double>("vy_gauss_mean");
        sim_conf.gauss_vystdd = cfg.get<double>("vy_gauss_stdd");
    } else {
        std::cerr << "Unknown initial y-velocity profile" << std::endl;
    }
    // Mass profile
    if( dmass=="uniform" ){
        sim_conf.distr_mass = uniform;
        sim_conf.unifo_mass_min = cfg.get<double>("m_uniform_min");
        sim_conf.unifo_mass_max = cfg.get<double>("m_uniform_max");
    } else if( dmass=="poisson" ){
        sim_conf.distr_mass = poisson;
        sim_conf.poiss_mass_mean = cfg.get<double>("m_poisson_mean");
    } else {
        std::cerr << "Unknown initial mass profile" << std::endl;
    }
    return 0;
}
