#include <cmath>
#include <mpi.h>
#include "iomanager.hpp"
#include "sroutines.hpp"
#include "nodes.hpp"
#include "config_file.hpp"

/** Initializes a I/O manager. The constructor sets the iteration counter of the
 * class, and calls an auxilary function to read the data provided in the
 * argument string.
 * @param[in] input_path Path to the configuration of the program. Note this
 * input path does not refer to any file. The function GenerateConfig will
 * append the expected extensions `.conf` and `.init` to it.
 */
IOManager::IOManager( std::string input_path )
    :write_iter(0)
{
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank==0 ) {
        std::cout << "CPU" << rank << " reads config..." << std::endl;
        GenerateConfig( input_path );
        std::cout << "CPU" << rank << " done reading!" << std::endl;
    }
    MPI_Barrier( MPI_COMM_WORLD );
    this->BroadcastConfig();
}

/** For a given Simulation object, writes the required data to output files.
 * These are formatted using the CSV format.
 * @param[in] sim Reference to the Simulation object containing the simulated
 * system. This obbject likely called the current method and passed itself by
 * reference.
 * @note In the current implementation, the user has no control about the output
 * (except the path it is written in). Future versions may allow selecting which
 * data is output, as well as the output frequency.
 */
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

/** Returns a constant reference to the current configuration.
 * Since we expect that the configuration object is initialized in the
 * constructor, no verification is performed before communicating it.
 */
const SConfig& IOManager::GetConfig() const {
    return conf_;
}

SError IOManager::DistributeTree( const Simulation& sim ) {
    int rank,size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    // Compute buffer sizes for particle transfer
    int bsize[2];
    // Implement uniform distribution of particles to start with
    bsize[0] = std::floor( (double)conf_.npart / (double)size );
    bsize[1] = bsize[0] + conf_.npart%size;
    if( rank==0 ) {
        std::cout << "CPU" << rank
            << " sends particles that other CPU should expect."
            << std::endl;
        Node *ptr( sim.tree_->GetNextLeaf( sim.tree_->GetRoot() ) );
        for( unsigned iskip(0); iskip<bsize[0]; iskip++ ) {
            std::cout << "  " << rank << ": Skipping " << *(ptr->GetParticle())
                << std::endl;
            ptr = sim.tree_->GetNextLeaf( ptr );
        }
        unsigned ip(1); // Starting index of receiving CPU
        while( ip<size ){
            unsigned send_size(bsize[ip==size-1]);
            Particle buf[send_size];
            unsigned pidx(0);
            while( ptr != NULL and pidx<send_size ){
                ///@todo Fix this absolute bullshit. We precisely implemented
                ///the operator= to remove the identifier.
                unsigned pid( ptr->GetParticle()->id );
                buf[pidx] = *(ptr->GetParticle());
                buf[pidx].id = pid;
                std::cout << "  " << rank << ": Buffering " << buf[pidx]
                    << std::endl;
                ptr = sim.tree_->GetNextLeaf( ptr );
                pidx++;
            }
            MPI_Request req;
            MPI_Isend( buf, send_size, MPI_Particle, ip, 0,
                    MPI_COMM_WORLD, &req );
            ip++;
        }
    } else {
        // File requests to get Particles
        unsigned recv_size(bsize[rank==size-1]);
        Particle buf[recv_size];
        MPI_Status s;
        MPI_Recv( buf, recv_size, MPI_Particle, 0, 0, MPI_COMM_WORLD, &s );
        std::cout << "CPU" << rank << " received " << bsize[rank==size-1] << " particles from CPU0." << std::endl;
        for( auto& part : buf ) {
            std::cout << "  " << rank << ": " << part << std::endl;
            Particle* newpart = new Particle(part);
            conf_.parts.push_back( newpart );
        }
    }
    // Delete removed particles from process zero
    MPI_Barrier( MPI_COMM_WORLD );
    if( rank==0 ) {
        std::cout << "Deleting memory from CPU0." << std::endl;
        Node* ptr( sim.tree_->GetNextLeaf( sim.tree_->GetRoot() ) );
        for( unsigned iskip(0); iskip<bsize[0]; iskip++ )
            ptr = sim.tree_->GetNextLeaf( ptr );
        while( ptr != NULL ) {
            std::cout << "  " << rank << ": Removing from memory " <<
                *(ptr->GetParticle()) << std::endl;
            std::vector<Particle*>::iterator it( conf_.parts.begin() );
            while( it != conf_.parts.end() ) {
                if( (*it)->id == ptr->GetParticle()->id ) {
                    delete *it;
                    conf_.parts.erase( it );
                    break;
                }
                it++;
            }
            ptr = sim.tree_->GetNextLeaf( ptr );
        }
    }
    return E_SUCCESS;
}

/** @brief Given an input path, creates a configuration corresponding to the
 * requested parameters.
 * @details This function is private and called by the constructor of the class
 * at initialization. It expects a path ending with the name of the simulation
 * to be run. Say this path reads `input/sim`, the function will look for the
 * files `input/sim.conf` and `input/sim.init` containing key-value parameters
 * and initial conditions respectively. This function dynamically allocates
 * memory for Particle objects that are read in the initial conditions file, but
 * note the class is not repsonsible for managing them: the SConfig structure is
 * responsible for destroying these when it reaches its end-of-life.
 * @param[in] file String indicating the path where the configuration and
 * initial conditions files can be found.
 */
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

SError IOManager::BroadcastConfig() {
    // Define buffers for (unsigned int) values of SConfig
    unsigned uv[2];
    double   dv[5];
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank==0 ) {
        uv[0] = conf_.npart;
        uv[1] = conf_.max_iter;
        dv[0] = conf_.dsize;
        dv[1] = conf_.dt;
        dv[2] = conf_.theta;
        dv[3] = conf_.max_wtime;
        dv[4] = conf_.extra_time;
        std::cout << "CPU" << rank << " broadcasts configuration." << std::endl;
    }
    MPI_Bcast( &uv, 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dv, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if( rank>0 ) {
        conf_.dsize = dv[0];
        conf_.npart = uv[0];
        conf_.theta = dv[2];
        conf_.max_iter = uv[1];
        conf_.max_wtime= dv[3];
        conf_.extra_time=dv[4];
        std::cout << "CPU" << rank << " received configuration." << std::endl;
    }
    return E_SUCCESS;
}

/** Wrapper function for opening a stream to an output file.
 * @param[in,out] ofstrm Output stream to a file.
 * @param[in] ofile Name of the file where the stream needs to point.
 */
SError IOManager::OpenStream( std::ofstream& ofstrm, const std::string& ofile ) const {
    ofstrm.open( ofile, std::ofstream::out );
    return E_SUCCESS;
}

/** Flushes an open stream and closes it.
 * @param[in] ofstrm Stream pointing to an output file.
 */
SError IOManager::CloseStream( std::ofstream& ofstrm ) const {
    ofstrm.flush();
    ofstrm.close();
    return E_SUCCESS;
}
