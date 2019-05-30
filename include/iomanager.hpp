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

class IOManager {
    public:
        IOManager( std::string input_path );
        SError WriteOutput( const Simulation& sim, unsigned it );
        const SConfig& GetConfig() const;
    private:
        SError GenerateConfig( const std::string& cfile, const std::string& ifile );
        SError OpenStream( std::ofstream& ofstrm, const std::string& ofile ) const;
        SError CloseStream( std::ofstream* ptr ) const;
        SConfig conf_;
};

#endif // SURMISE_IOMANAGER_HPP_
