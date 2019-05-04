#include <unistd.h>
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

int ReadInputFile( SConfig& sim_conf, char* filename ) {
    sim_conf.yo = 10;
    sim_conf.str = "Punchline bb!";
    return 0;
}
