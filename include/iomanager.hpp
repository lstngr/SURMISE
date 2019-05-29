/**
 * @file IOManager
 */

#ifndef SURMISE_IOMANAGER_HPP_
#define SURMISE_IOMANAGER_HPP_

#include <fstream>
#include <string>
#include "serrors.hpp"
#include "stypes.hpp"

class IOManager {
    public:
        IOManager( std::string filename );
        IOManager( );
        bool IsReady() const;
    protected:
    private:
        SError OpenStream();
        SError CloseStream();
        std::ostream* outstream_;
};

#endif // SURMISE_IOMANAGER_HPP_
