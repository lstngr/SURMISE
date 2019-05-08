/** \file sinput.hpp
 * @brief Contains routines managing reading the input file containing the
 * simulation parameters. Most of the relevant code reading the file is
 * available in config_file.hpp.
 */
#ifndef SURMISE_SINPUT_HPP_
#define SURMISE_SINPUT_HPP_
#include "stypes.hpp"

/** @brief Reads the arguments provided with the call to the binary.
 * @param[in] error_code An error code that will be non-zero if the function
 * fails.
 * @return Structure containing the relevant parameters.
 */
SConfig ReadInputs( int argc, char** argv, int& error_code );

/** @brief Fills the simulation configuration structure with input parameters
 * read from a parameter file.
 * @param[in,out] sim_conf Simulation configuration that will be overwritten with
 * the values read from filename.
 * @param[in] filename Name of the file containing the input parameters.
 * @param[in] filename_init Name of the file containing the initial particle
 * distributions
 * @return Error code. Non-zero if the routine failed.
 */
int ReadInputFile( SConfig& sim_conf, const std::string& filename, const std::string& filename_init );

#endif //SURMISE_SINPUT_HPP_

