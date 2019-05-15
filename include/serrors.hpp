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
    /// Indicates a routine was executed as planned and exited without error.
    //Conventionally, we give it the value zero.
    E_SUCCESS = 0
};

#endif // SURMISE_SERRORS_HPP_
