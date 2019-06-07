#ifndef SURMISE_TIMER_HPP_
#define SURMISE_TIMER_HPP_

#include <iostream>
#include "serrors.hpp"
#include "iomanager.hpp"

/** Enumeration assigning indices to timer tasks.
 * @warning Relies on the last element to get the number of elements defined in
 * the enumeration. Behavior may vary across compilers.
 */
enum STimerNames {
    T_ITER=0,
    T_DISTR,
    T_FORCES,
    T_EVOL,
    T_LEAFSYNC,
    T_TREEUPDATE,
    T_OUTPUT,
    T_COUNT //!< Count of the number of elements in STimerNames (may vary across compilers)
};

/** A class measuring the almighty time. It relies exclusively on MPI_Wtime.
 */
class Timer {
    public:
        Timer( unsigned count=T_COUNT );
        ~Timer();
        SError StartTimer( unsigned idx );
        SError SetTimer( unsigned idx, const double& val );
        SError StopTimer(  unsigned idx );
        bool IsActive( unsigned idx ) const;
        unsigned GetNumber() const;
        double GetTime( unsigned idx ) const;
    private:
        /** Integer indicating the number of timers in the class.
         */
        unsigned counters_;
        /** Array of times. Times computed between a timer start and stop calls
         * are stored here.
         */
        double* timers_;
        /** Array of start times. These are temporary, are set by a start call,
         * and are erased by a stop call.
         */
        double* starts_;
};
std::ostream& operator<<( std::ostream& os, const Timer& t );

#endif // SURMISE_TIMER_HPP_
