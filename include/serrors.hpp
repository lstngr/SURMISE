/** @file serrors.hpp Defines error code for the SURMISE code.
 * Details
 */

#ifndef SURMISE_SERRORS_HPP_
#define SURMISE_SERRORS_HPP_

/** @brief Contains all error codes within the SURMISE project.
 * @details This enumeration is used to manage errors within the SURMISE code.
 * It associates a keyword linked to the encountered error and binds it to an
 * integer value, which can easily be return by most methods.
 * All errors are prefixed by `E_`.
 */
enum SError {
    /** Indicates a routine was executed as planned and exited without error.
    * Conventionally, we give it the value zero.*/
    E_SUCCESS = 0,
    /** Returned to indicate a general failure. Execution is aborted if
     * encountered.*/
    E_FAILURE = 1,
    /** Returned when a given task exceeds the time that was allocated for it to
     * complete. For example, if the main Simulation::Run routine still runs
     * when the maximum allocated amount of time in the configuration is
     * reached.*/
    E_TIMEOUT = 2
};

#endif // SURMISE_SERRORS_HPP_
