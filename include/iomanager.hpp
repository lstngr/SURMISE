/**
 * @file IOManager
 */

#ifndef SURMISE_IOMANAGER_HPP_
#define SURMISE_IOMANAGER_HPP_

#include <fstream>
#include <string>
#include "serrors.hpp"
#include "stypes.hpp"

class Simulation;

/** @brief Manager for Input and Output Routines of the SURMISE code.
 * @details This class is responsible for reading the program's intialization
 * files at startup, as well as exporting the requested results to output files.
 * This class also initializes the configuration structure, of type SConfig,
 * which contains the dynamically allocated Particle structures, based on the
 * initial condition file provided.
 */
class IOManager {
    public:
        IOManager( std::string input_path );
        SError WriteOutput( const Simulation& sim );
        const SConfig& GetConfig() const;
        SError DistributeParticles( Simulation& sim );
        SError DistributeTree( const Simulation& sim );
        SError SyncLeafs( Simulation& sim );
    private:
        SError GenerateConfig( const std::string& file );
        SError BroadcastConfig();
        SError OpenStream( std::ofstream& ofstrm, const std::string& ofile ) const;
        SError CloseStream( std::ofstream& ofstrm ) const;
        SConfig conf_; //!< Object containing the parameters of the simulation.
        unsigned write_iter; //!< Counter labeling how many writes occured.
                             //!< Useful for naming output files.
};

#endif // SURMISE_IOMANAGER_HPP_
