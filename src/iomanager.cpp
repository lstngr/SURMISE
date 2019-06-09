#include <cmath>
#include <iomanip>
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
    int rank,size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    std::ofstream particles, tree, timers;
    OpenStream( particles, conf_.opath + ".leafs." + std::to_string(write_iter) );
    OpenStream( tree, conf_.opath + ".tree." + std::to_string(write_iter) );
    if( sim.timer_ != NULL )
        OpenStream( timers, conf_.opath + ".timers." + std::to_string(write_iter) );
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
    // TIMER PRINTING
    if( sim.timer_ != NULL ) {
        double *sbuf(NULL), *rbuf(NULL);
        sbuf = new double[T_COUNT];
        // Get all timers on CPU0
        rbuf = new double[size*T_COUNT];
        for( unsigned i_timer(0); i_timer<T_COUNT; i_timer++ )
            sbuf[i_timer] = sim.timer_->GetTime(i_timer);
        MPI_Gather( sbuf, T_COUNT, MPI_DOUBLE, rbuf, T_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD );
        delete[] sbuf;
        for( int ip(0); ip<size; ip++ ) {
            timers << ip << std::setprecision(8);
            for( unsigned it(0); it<T_COUNT; it++ ) {
                timers << "," << rbuf[ip*T_COUNT+it];
            }
            timers << std::endl;
        }
        delete[] rbuf;
    }
    CloseStream( particles );
    CloseStream( tree );
    if( sim.timer_ != NULL )
        CloseStream( timers );
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

/** @brief Sends a copy of Particles stored in the master process and builds an
 * update list for each thread.
 * @details The routine is called by all (participating) threads (usually, all
 * running instances). It assumes the run just began, and that only the
 * configuration object of the master process is filled with particles.
 * First, particles are being broadcast from the master process, and indexes are
 * manipulated so that they're preserved (remember the assignement operator
 * resets it otherwise).
 * Then, each thread builds a tree local to its process.
 * Finally, the number of particles handled by each processor is computed, and
 * the leaf level of the tree is browsed with QuadTree::GetNextLeaf. Booleans
 * are set to indicate a process whether it manages a given particle or not. The
 * ordering of the Z-curve should pack particles close to another and assign
 * them to a similar process.
 *
 * @param[in] sim Reference to the simulation object to modify. It is used for
 * Simulation related tasks, such as tree building and browsing.
 * @returns An error code.
 */
SError IOManager::DistributeParticles( Simulation& sim ) {
    int rank,size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    Particle *buf = new Particle[conf_.npart];
    if( rank==0 ) {
        for( unsigned ip(0); ip<conf_.npart; ip++ ) {
            buf[ip] = *(conf_.parts[ip]);
            buf[ip].id = conf_.parts[ip]->id;
        }
    }
    MPI_Bcast( buf, conf_.npart, MPI_Particle, 0, MPI_COMM_WORLD );
    if( rank>0 ) {
        for( unsigned i(0); i<conf_.npart; i++ ) {
            const Particle& p( buf[i] );
            conf_.parts.push_back( new Particle(p) );
            conf_.parts.back()->id = p.id;
        }
    }
    delete[] buf;
    sim.BuildTree();
    unsigned upsize[2];
    upsize[0] = std::floor( (double)conf_.npart / (double)size );
    upsize[1] = conf_.npart - (size-1)*upsize[0];
    sim.update_list_ = new bool[conf_.npart];
    Node* ptr( sim.tree_->GetNextLeaf( sim.tree_->GetRoot() ) );
    unsigned ic(0);
    while(ptr!=NULL) {
        unsigned ip(ptr->GetParticle()->id);
        sim.update_list_[ip]=false;
        if(rank*upsize[0] <= ic and (rank+1)*upsize[0] > ic and rank<size-1) {
            sim.update_list_[ip] = true;
        } else if( ic>=rank*upsize[0] and rank==size-1 ) {
            sim.update_list_[ip] = true;
        }
        ptr = sim.tree_->GetNextLeaf( ptr );
        ic++;
    }
    return E_SUCCESS;
}

/** @brief Ensures all processors share the most up-to-date information about
 * every particle in the system.
 * @details Given a simulation object, the threads first communicate how many
 * particles they each manage to one another. A buffer is then prepared, in
 * which particles managed by the running thread are stored (with their index
 * preserved). A nonblocking send request is then filed through the MPI
 * interface.
 * Then, each thread enter a loop where they wait for each processor to
 * communicate their particles. The receiving operation is blocking. Received
 * particles are stored in the configuration where they override their _older_
 * counterparts.
 * The routine finally frees all allocated buffers in memory and returns.
 *
 * @warning This routine assumes all particles are available on all processors,
 * and that they are ordered by increasing index in the configuration object. It
 * then uses the particle's index to avoid searching through local particles
 * which is the one to replace, but could lead to segmentation faults or logical
 * errors easily if other parts of the code change.
 *
 * @param[in] sim Reference to the simulation object in which the particle
 * update lists are stored.
 * @returns An error code.
 */
SError IOManager::SyncLeafs( Simulation& sim ) {
    int rank,size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    unsigned bsize(0), *bsizes, ib(0);
    Particle *recv_buf, *send_buf;
    bsizes = new unsigned[size];
    for(unsigned i(0); i<conf_.npart; i++){
        if( sim.update_list_[i] )
            bsize++;
    }
    MPI_Allgather( &bsize, 1, MPI_UNSIGNED, bsizes, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
    send_buf = new Particle[bsize];
    for(unsigned i(0); i<conf_.npart; i++){
        if( sim.update_list_[i] ){
            send_buf[ib] = *(conf_.parts[i]);
            send_buf[ib].id = conf_.parts[i]->id;
            ib++;
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Request *send_req = new MPI_Request[size];
    for( int i(0); i<size; i++ ){
        if( i!=rank ){
            MPI_Isend( send_buf, bsizes[rank], MPI_Particle, i, 0,
                    MPI_COMM_WORLD, &send_req[i] );
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );
    for( int i(0); i<size; i++ ){
        if( i!=rank ){
            recv_buf = new Particle[bsizes[i]];
            MPI_Status stat;
            MPI_Recv( recv_buf, bsizes[i], MPI_Particle, i, 0, MPI_COMM_WORLD,
                    &stat );
            for( unsigned is(0); is<bsizes[i]; is++ ) {
                unsigned identifier( recv_buf[is].id );
                *(conf_.parts[identifier]) = recv_buf[is];
                conf_.parts[identifier]->id = identifier;
                // for( const auto& pp : conf_.parts ) {
                //     if( pp->id == identifier ) {
                //         *pp = recv_buf[is];
                //         pp->id = identifier;
                //         break;
                //     }
                // }
            }
            delete[] recv_buf;
        }
    }
    delete[] bsizes;
    delete[] send_buf;
    delete[] send_req;
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
    conf_.npart      = cfg.get<unsigned>("npart");
    conf_.dt         = cfg.get<double>("tevol_dt");
    conf_.theta      = cfg.get<double>("theta");
    conf_.max_iter   = cfg.get<unsigned>("max_iter");
    conf_.max_wtime  = cfg.get<double>("walltime");
    conf_.extra_time = cfg.get<double>("extratime");
    conf_.write_freq = cfg.get<unsigned>("write_freq");
    if( not conf_.write_freq )
        conf_.write_freq = 1;
    std::cout << "Read write_freq=" << conf_.write_freq << std::endl;
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

/** Distributes all numeric parameters of the configuration object to available
 * threads.
 */
SError IOManager::BroadcastConfig() {
    // Define buffers for (unsigned int) values of SConfig
    unsigned uv[3];
    double   dv[5];
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( rank==0 ) {
        uv[0] = conf_.npart;
        uv[1] = conf_.max_iter;
        uv[2] = conf_.write_freq;
        dv[0] = conf_.dsize;
        dv[1] = conf_.dt;
        dv[2] = conf_.theta;
        dv[3] = conf_.max_wtime;
        dv[4] = conf_.extra_time;
        std::cout << "CPU" << rank << " broadcasts configuration." << std::endl;
    }
    MPI_Bcast( &uv, 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dv, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    if( rank>0 ) {
        conf_.dsize = dv[0];
        conf_.npart = uv[0];
        conf_.dt = dv[1];
        conf_.theta = dv[2];
        conf_.max_iter = uv[1];
        conf_.max_wtime= dv[3];
        conf_.extra_time=dv[4];
        conf_.write_freq=uv[2];
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
