#include "iomanager.hpp"
#include "config_file.hpp"

IOManager::IOManager( std::string input_path )
{
    std::string conf_path( input_path + ".conf" );
    std::string init_path( input_path + ".init" );
    GenerateConfig( conf_path, init_path );
}

SError IOManager::WriteOutput( const Simulation& sim, unsigned it ) {
    std::ofstream* ofstrm_conf(NULL);
    OpenStream( *ofstrm_conf, conf_.opath + ".conf." + std::to_string(it) );
    // Do something here
    CloseStream( ofstrm_conf );
    return E_SUCCESS;
}

const SConfig& IOManager::GetConfig() const {
    return conf_;
}

SError IOManager::GenerateConfig( const std::string& cfile, const std::string& ifile ) {
    ConfigFile cfg(cfile);
    // Load basic parameters
    conf_.dsize      = cfg.get<double>("domsize");
    conf_.npart      = cfg.get<int>("npart");
    conf_.dt         = cfg.get<double>("tevol_dt");
    conf_.epsilon    = cfg.get<double>("epsilon");
    conf_.max_iter   = cfg.get<int>("max_iter");
    conf_.max_wtime  = cfg.get<double>("walltime");
    conf_.extra_time = cfg.get<double>("extratime");
    conf_.opath      = cfg.get<std::string>("output_path");
    // Load particles
    // Particle array is dynamically allocated and should match input file size.
    conf_.parts.reserve(conf_.npart);
    std::ifstream ifsp( ifile );
    for( size_t ip(0); ip < conf_.npart; ip++  ) {
        conf_.parts.push_back( new Particle );
        conf_.parts.back()->id = ip;
        ifsp >> conf_.parts.back()->pos[0]
             >> conf_.parts.back()->pos[1]
             >> conf_.parts.back()->vel[0]
             >> conf_.parts.back()->vel[1]
             >> conf_.parts.back()->mass;
        conf_.parts.back()->frc[0] = 0.0;
        conf_.parts.back()->frc[1] = 0.0;
    }
    ifsp.close();
    return E_SUCCESS;
}

SError IOManager::OpenStream( std::ofstream& ofstrm, const std::string& ofile ) const {
    return E_SUCCESS;
}
