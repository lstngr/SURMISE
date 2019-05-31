#include "iomanager.hpp"
#include "sroutines.hpp"
#include "nodes.hpp"
#include "config_file.hpp"

IOManager::IOManager( std::string input_path )
    :write_iter(0)
{
    GenerateConfig( input_path );
}

SError IOManager::WriteOutput( const Simulation& sim ) {
    std::ofstream particles;
    OpenStream( particles, conf_.opath + ".leafs." + std::to_string(write_iter) );
    // Get Root of the tree
    Node* ptr(sim.tree_->GetRoot());
    ptr = sim.tree_->GetNextLeaf( ptr );
    // Iterate over tree (GetNextLeaf) until found all paptricles.
    while( ptr!=NULL ) {
        particles << *ptr << std::endl;
        ptr = sim.tree_->GetNextLeaf( ptr );
    }
    CloseStream( particles );
    write_iter++;
    return E_SUCCESS;
}

const SConfig& IOManager::GetConfig() const {
    return conf_;
}

SError IOManager::GenerateConfig( const std::string& file ) {
    std::string cfile( file+".conf" );
    std::string ifile( file+".init" );
    std::string oname( file );
    oname = oname.substr( oname.find_last_of("/\\")+1, oname.length());
    ConfigFile cfg(cfile);
    // Load basic parameters
    conf_.dsize      = cfg.get<double>("domsize");
    conf_.npart      = cfg.get<int>("npart");
    conf_.dt         = cfg.get<double>("tevol_dt");
    conf_.theta      = cfg.get<double>("theta");
    conf_.max_iter   = cfg.get<int>("max_iter");
    conf_.max_wtime  = cfg.get<double>("walltime");
    conf_.extra_time = cfg.get<double>("extratime");
    conf_.ipath      = file;
    conf_.opath      = cfg.get<std::string>("output_path") + oname;
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
    ofstrm.open( ofile, std::ofstream::out );
    return E_SUCCESS;
}

SError IOManager::CloseStream( std::ofstream& ofstrm ) const {
    ofstrm.flush();
    ofstrm.close();
    return E_SUCCESS;
}
