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
    std::ofstream particles, tree;
    OpenStream( particles, conf_.opath + ".leafs." + std::to_string(write_iter) );
    OpenStream( tree, conf_.opath + ".tree." + std::to_string(write_iter) );
    // Get Root of the tree
    Node* pptr(sim.tree_->GetRoot());
    pptr = sim.tree_->GetNextLeaf( pptr );
    Node* tptr(sim.tree_->GetRoot());
    // PARTICLE PRINTING
    // Iterate over tree (GetNextLeaf) until found all papptricles.
    while( pptr!=NULL ) {
        particles << pptr->com_->id << "," << pptr->com_->mass << ","
            << pptr->com_->pos[0] << "," << pptr->com_->pos[1] << ","
            << pptr->com_->vel[0] << "," << pptr->com_->vel[1] << ","
            << pptr->com_->frc[0] << "," << pptr->com_->frc[1]
            << std::endl;
        pptr = sim.tree_->GetNextLeaf( pptr );
    }
    // TREE PRINTING
    tree << tptr << "," << tptr->children_[0] << ","
        << tptr->children_[1] << ","
        << tptr->children_[2] << ","
        << tptr->children_[3] << std::endl;
    bool print_tree(true);
    while(print_tree) {
        // Goes as deep as it can, prints all nodes when stepping down.
        if( not tptr->IsLeaf() ) {
            for( unsigned ic(0); ic<4; ic++ ) {
                if( tptr->GetChild(ic) != NULL ) {
                    tptr = tptr->GetChild(ic);
                    tree << tptr << "," << tptr->children_[0] << ","
                         << tptr->children_[1] << ","
                         << tptr->children_[2] << ","
                         << tptr->children_[3] << std::endl;
                    break;
                }
            }
        } else {
            // We reached a leaf. Need to go up until we can move again.
            bool go_up(true);
            while(go_up){
                // If we're at the root, cannot climb
                if(tptr->IsRoot()){
                    print_tree=false;
                    break;
                }
                // Else get the parent.
                // Figure out our index from parent
                unsigned icur(tptr->GetParent()->GetQuadrant(tptr));
                tptr = tptr->GetParent();
                for( unsigned ic(icur+1); ic<4; ic++ ) {
                    if( tptr->GetChild(ic) != NULL ){
                        // Found a child from parent
                        tptr = tptr->GetChild(ic);
                        tree << tptr << "," << tptr->children_[0] << ","
                            << tptr->children_[1] << ","
                            << tptr->children_[2] << ","
                            << tptr->children_[3] << std::endl;
                        go_up = false;
                        break;
                    }
                }
                // Here, we're about to redo the while.
                // If go_up = false, we found a child to go to.
                // Else, the routine restarts and goes up one other level.
            }
        }
    }
    CloseStream( particles );
    CloseStream( tree );
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
